% --- generate_test_dataset.m ---
%
% Generates a TEST dataset with UNKNOWN SEEDS (different PN initial states).
%
% Output: data/synthetic/classification_dataset_stft_awgn_test.mat

clear; clc; close all;

addpath(genpath('src'));
fprintf('Added src/ and subfolders to the MATLAB path.\n');

% --- Simulation + dataset configuration ---
snr_levels_db = -10:5:10; % Reduced SNR granularity for faster testing (-10, -5, 0, 5, 10)
num_signals_per_config = 1; % Reduced for faster testing

% Gold code setup (degree-5 preferred pair)
pnPoly1 = [5 2 0];
pnPoly2 = [5 4 3 2 0];

% NEW Initial states (NOT in training set)
hopset_initial_states_test = { ...
    [1 1 1 1 0], [0 1 1 1 0], [0 1 0 1 0], [1 0 1 0 0], [0 0 1 0 1] ...
    };

% STFT parameters
stft_params.windowLength = 128;
stft_params.overlapLength = 96; % 75% overlap
stft_params.nfft = 1024;

% Chunk write size
chunkSize = 4096;

% --- Define FHSS simulation parameters ---
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

% FSK modulator
modulator = comm.FSKModulator(simParams.M, simParams.freqSeparation, ...
    'SamplesPerSymbol', simParams.samplesPerSymbol, ...
    'SymbolRate', simParams.symbolRate);

% --- MAT-file initialization ---
output_filename = 'data/synthetic/classification_dataset_stft_awgn_test.mat';
if exist(output_filename, 'file'), delete(output_filename); end
m = matfile(output_filename, 'Writable', true);

total_signals = numel(snr_levels_db) * numel(hopset_initial_states_test) * num_signals_per_config;
total_examples = total_signals * simParams.numHops;

% Pre-allocate variables
m.X_data = cell(total_examples, 1);
m.Y_data = zeros(total_examples, 1, 'uint8');
m.numHopsPerSignal = simParams.numHops;
m.stft_params = stft_params;

signal_counter = 0;
example_write_counter = 0;
chunk_idx = 0;
X_chunk = cell(chunkSize, 1);
Y_chunk = zeros(chunkSize, 1, 'uint8');

fprintf('Starting TEST dataset generation...\n');
fprintf('Total signals: %d | Total hops: %d\n', total_signals, total_examples);

for i = 1:numel(snr_levels_db)
    snr_db = snr_levels_db(i);
    fprintf('SNR = %d dB\n', snr_db);

    for j = 1:numel(hopset_initial_states_test)
        simParams.pnInitial2 = hopset_initial_states_test{j};

        for k = 1:num_signals_per_config
            % --- Generate hop sequence ---
            goldCodeSequence = generateGoldCodeSequence(simParams);
            reshapedBits = reshape(goldCodeSequence, simParams.k, []);
            transposedBits = reshapedBits';
            hopIndices_col = bi2de(transposedBits, 'left-msb');
            hopIndices = hopIndices_col';
            hopFrequencies = simParams.hopset(hopIndices + 1);

            % --- Generate modulated data ---
            messageBits = randi([0 1], simParams.numBits, 1);
            modulatedData = modulator(messageBits);

            % --- Per-hop processing ---
            for hop_idx = 1:simParams.numHops
                chunk_idx = chunk_idx + 1;

                startIndex = (hop_idx-1) * simParams.samplesPerHop + 1;
                stopIndex = hop_idx * simParams.samplesPerHop;
                dataSegment = modulatedData(startIndex:stopIndex);
                t = (0:length(dataSegment)-1)' / simParams.fs;
                hoppingCarrier = exp(1j * 2*pi * hopFrequencies(hop_idx) * t);
                singleHopSignal = dataSegment .* hoppingCarrier;

                receivedHopSignal = addAWGN(singleHopSignal, snr_db);
                [S, ~, ~] = spectrogram(receivedHopSignal, hamming(stft_params.windowLength), ...
                    stft_params.overlapLength, stft_params.nfft, simParams.fs, "centered");

                X_chunk{chunk_idx} = single(20*log10(abs(S) + eps));
                Y_chunk(chunk_idx) = uint8(hopIndices(hop_idx));

                if chunk_idx == chunkSize
                    start_write_idx = example_write_counter + 1;
                    end_write_idx = example_write_counter + chunkSize;
                    m.X_data(start_write_idx:end_write_idx, 1) = X_chunk;
                    m.Y_data(start_write_idx:end_write_idx, 1) = Y_chunk;
                    example_write_counter = end_write_idx;
                    chunk_idx = 0;
                end
            end

            signal_counter = signal_counter + 1;
        end
    end
end

% Write any remaining chunk
if chunk_idx > 0
    start_write_idx = example_write_counter + 1;
    end_write_idx = example_write_counter + chunk_idx;
    m.X_data(start_write_idx:end_write_idx, 1) = X_chunk(1:chunk_idx);
    m.Y_data(start_write_idx:end_write_idx, 1) = Y_chunk(1:chunk_idx);
    example_write_counter = end_write_idx;
end

% --- Calculate Global Normalization Stats ---
% We iterate through the file to be memory-safe, or we could have tracked it above.
% Since we just wrote it, let's track it in the loop above?
% Actually, re-reading or tracking during loop is better.
% Let's modify the loop structure in a future step if needed, but for now,
% since we already wrote the code, let's just add a quick pass or rely on the python script
% calculating it dynamically as it does currently.

% HOWEVER, the user explicitly asked to "normalize correctly".
% The best way is to Calculate it during generation.

% Let's rewrite the Min/Max section to be computed *during* the generation loop
% I will patch the loop range in the previous block.
% But since I can't easily jump back and edit the loop variable initialization in this single replace call without replacing the whole file,
% I will append a calculation block here that reads the data in chunks (since it is already saved).

fprintf('Calculating global min/max for normalization...\n');
global_min = Inf;
global_max = -Inf;

% Re-read in chunks to find min/max
num_chunks = ceil(example_write_counter / chunkSize);
for i = 1:num_chunks
    s_idx = (i-1)*chunkSize + 1;
    e_idx = min(i*chunkSize, example_write_counter);

    % partial load
    data_cells = m.X_data(s_idx:e_idx, 1);

    for k = 1:numel(data_cells)
        val = data_cells{k};
        if ~isempty(val)
            max_val = max(val, [], 'all');
            min_val = min(val, [], 'all');
            if max_val > global_max, global_max = max_val; end
            if min_val < global_min, global_min = min_val; end
        end
    end
end

m.global_min = global_min;
m.global_max = global_max;
fprintf('Global Min: %.4f | Global Max: %.4f\n', global_min, global_max);

fprintf('TEST Dataset generation complete. Saved %d examples to %s\n', example_write_counter, output_filename);
