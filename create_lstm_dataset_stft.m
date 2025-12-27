% --- create_lstm_dataset_stft.m ---
%
% Generates a memory-efficient, batched dataset of FHSS signals.
%
% This version uses a matfile object and batching to handle datasets that
% are too large to fit in RAM. It writes data to disk in chunks to
% balance memory usage and performance.
%
% Output file: 'data/synthetic/classification_dataset_stft_awgn.mat'
%   - X_data: A cell array where each cell is a [Frequency x Time] matrix
%             representing the log-magnitude STFT of a *single* hop.
%   - Y_data: A vector where each element is the hop index (label)
%             corresponding to the STFT in X_data.

clear; clc; close all;

% --- 0. Setup Project Path ---
addpath(genpath('src'));
fprintf('Added src/ and all subfolders to the MATLAB path.\n');

% --- 1. Configuration ---

% SNR levels in dB
snr_levels_db = -10:2:10;

% Number of signals to generate per hopset and SNR level
num_signals_per_config = 3;

% Gold Code "preferred pair" polynomials of degree 5
pnPoly1 = [5 2 0];
pnPoly2 = [5 4 3 2 0];

% Define different hopsets
hopset_initial_states_g2 = { ...
    [0 0 0 0 1], [0 0 0 1 0], [0 0 1 1 1], [0 1 0 0 1], [0 1 1 0 0], ...
    [1 0 0 0 1], [1 0 1 1 0], [1 1 0 0 1], [1 1 1 0 1], [1 1 1 1 1] ...
    };
num_hopsets = length(hopset_initial_states_g2);

% STFT Parameters
stft_params = struct();
stft_params.windowLength = 128;   %%%
stft_params.overlapLength = 96; % 75% overlap
stft_params.nfft = 1024;

% Batching configuration
chunkSize = 4096; % Number of examples to batch in memory before writing

% --- 2. Define Base Simulation Parameters ---
simParams = struct();
simParams.fs = 10e6;
simParams.numBits = 512;
simParams.M = 2;
simParams.symbolRate = simParams.fs * 0.01;
simParams.freqSeparation = simParams.symbolRate;
simParams.samplesPerSymbol = round(simParams.fs / simParams.symbolRate);
simParams.k = 3;
simParams.numHops = 256;
simParams.bitsPerFrame = simParams.numHops * simParams.k;
numChannels = 2^simParams.k;
spacing = 2*simParams.freqSeparation;
baseFreq = 2e6;
simParams.hopset = (0:numChannels-1) * spacing + baseFreq;
simParams.bitsPerSymbol = log2(simParams.M);
simParams.numSymbols = simParams.numBits / simParams.bitsPerSymbol;
simParams.numModulatedSamples = simParams.numSymbols * simParams.samplesPerSymbol;
simParams.samplesPerHop = floor(simParams.numModulatedSamples / simParams.numHops);
simParams.pnPoly1 = pnPoly1;
simParams.pnInitial1 = [0 0 0 0 1];
simParams.pnPoly2 = pnPoly2;

% FSK Modulator (defined once for efficiency)
modulator = comm.FSKModulator(simParams.M, simParams.freqSeparation, ...
    'SamplesPerSymbol', simParams.samplesPerSymbol, ...
    'SymbolRate', simParams.symbolRate);

% --- 3. Dataset Initialization & matfile setup ---
output_filename = 'data/synthetic/classification_dataset_stft_awgn.mat';

if exist(output_filename, 'file'), delete(output_filename); end
m = matfile(output_filename, 'Writable', true);

total_signals = length(snr_levels_db) * num_hopsets * num_signals_per_config;
total_examples = total_signals * simParams.numHops;

% Pre-allocate variables in the MAT-file.
m.X_data = cell(total_examples, 1);
m.Y_data = zeros(total_examples, 1);

% Initialize counters and in-memory chunks
signal_counter = 0;
example_write_counter = 0;
chunk_idx = 0;
X_chunk = cell(chunkSize, 1);
Y_chunk = zeros(chunkSize, 1);

global_min = inf;
global_max = -inf;

% --- 4. Data Generation Loop ---
fprintf('Starting dataset generation...\n');
fprintf('Total examples to generate: %d\n', total_examples);

for i = 1:length(snr_levels_db)
    snr_db = snr_levels_db(i);
    fprintf('Processing SNR = %d dB\n', snr_db);

    for j = 1:num_hopsets
        simParams.pnInitial2 = hopset_initial_states_g2{j};

        for k = 1:num_signals_per_config
            % --- A. Generate Hop Sequence ---
            goldCodeSequence = generateGoldCodeSequence(simParams);
            reshapedBits = reshape(goldCodeSequence, simParams.k, []);
            % Break down the complex conversion to avoid syntax errors
            transposedBits = reshapedBits';
            hopIndices_col = bi2de(transposedBits, 'left-msb');
            hopIndices = hopIndices_col'; % Ensure hopIndices is a row vector
            hopFrequencies = simParams.hopset(hopIndices + 1);

            % --- B. Generate Modulated Data ---
            messageBits = randi([0 1], simParams.numBits, 1);
            modulatedData = modulator(messageBits);

            % --- C. Process and Store Data for Each Hop ---
            for hop_idx = 1:simParams.numHops
                chunk_idx = chunk_idx + 1;

                % 1. Generate signal for a single hop
                startIndex = (hop_idx-1)*simParams.samplesPerHop + 1;
                stopIndex = hop_idx*simParams.samplesPerHop;
                dataSegment = modulatedData(startIndex:stopIndex);
                t = (0:length(dataSegment)-1)' / simParams.fs;
                hoppingCarrier = exp(1j*2*pi * hopFrequencies(hop_idx) * t);
                singleHopSignal = dataSegment .* hoppingCarrier;

                % 2. Add noise and compute STFT
                receivedHopSignal = addAWGN(singleHopSignal, snr_db);
                [S, ~, ~] = spectrogram(receivedHopSignal, hamming(stft_params.windowLength), ...
                    stft_params.overlapLength, stft_params.nfft, simParams.fs, "centered");

                % 3. Store in chunk
                X_chunk{chunk_idx} = 20*log10(abs(S));
                Y_chunk(chunk_idx) = hopIndices(hop_idx);

                % 4. Write chunk to disk if full
                if chunk_idx == chunkSize
                    start_write_idx = example_write_counter + 1;
                    end_write_idx = example_write_counter + chunkSize;

                    m.X_data(start_write_idx:end_write_idx, 1) = X_chunk;
                    m.Y_data(start_write_idx:end_write_idx, 1) = Y_chunk;

                    % Update global stats
                    % Convert cell array of matrices to a large tensor/matrix for fast min/max
                    % X_chunk is {128xN, ...}. cat(3, ...) makes it 128xNx4096.
                    chunk_data_mat = cat(3, X_chunk{:});
                    current_min = min(chunk_data_mat(:));
                    current_max = max(chunk_data_mat(:));

                    global_min = min(global_min, current_min);
                    global_max = max(global_max, current_max);

                    example_write_counter = end_write_idx;
                    chunk_idx = 0; % Reset chunk index
                end
            end

            signal_counter = signal_counter + 1;
            if mod(signal_counter, 20) == 0
                fprintf('    Generated %d/%d signals. Total examples written: %d/%d\n', signal_counter, total_signals, example_write_counter, total_examples);
            end
        end
    end
end

% --- 5. Write Final Partial Chunk ---
if chunk_idx > 0
    start_write_idx = example_write_counter + 1;
    end_write_idx = example_write_counter + chunk_idx;
    fprintf('  Writing final chunk to disk (%d to %d) ...', start_write_idx, end_write_idx);

    m.X_data(start_write_idx:end_write_idx, 1) = X_chunk(1:chunk_idx);
    m.Y_data(start_write_idx:end_write_idx, 1) = Y_chunk(1:chunk_idx);

    % Update global stats for the final chunk
    final_chunk_data = cat(3, X_chunk{1:chunk_idx}); % Concatenate to 3D array
    current_min = min(final_chunk_data(:));
    current_max = max(final_chunk_data(:));
    global_min = min(global_min, current_min);
    global_max = max(global_max, current_max);

    example_write_counter = end_write_idx;
end

% --- 6. Save Global Stats ---
m.global_min = global_min;
m.global_max = global_max;
fprintf('\nGlobal Stats: Min = %.4f, Max = %.4f\n', global_min, global_max);


% --- 6. Finalization ---
fprintf('Dataset generation complete. Saved %d total examples to %s.\n', ...
    example_write_counter, output_filename);