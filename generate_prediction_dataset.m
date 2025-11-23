% --- generate_prediction_dataset.m ---
%
% Builds the base per-hop STFT dataset for FHSS signals and stores it in a
% MAT-file using a matfile handle for chunked, low-memory writes.
%
% Output: data/synthetic/classification_dataset_stft_awgn.mat
%   - X_data: cell array, each cell is [freq x time] log-magnitude STFT for a hop
%   - Y_data: numeric hop index label for each hop
%   - numHopsPerSignal: number of hops in each signal (for downstream scripts)
%   - stft_params: struct with STFT config (windowLength, overlapLength, nfft)

clear; clc; close all;

addpath(genpath('src'));
fprintf('Added src/ and subfolders to the MATLAB path.\n');

% --- Simulation + dataset configuration ---
snr_levels_db = -10:2:10;
num_signals_per_config = 3;

% Gold code setup (degree-5 preferred pair)
pnPoly1 = [5 2 0];
pnPoly2 = [5 4 3 2 0];

hopset_initial_states_g2 = { ...
    [0 0 0 0 1], [0 0 0 1 0], [0 0 1 1 1], [0 1 0 0 1], [0 1 1 0 0], ...
    [1 0 0 0 1], [1 0 1 1 0], [1 1 0 0 1], [1 1 1 0 1], [1 1 1 1 1] ...
};

% STFT parameters
stft_params.windowLength = 128;
stft_params.overlapLength = 96; % 75% overlap
stft_params.nfft = 1024;

% Chunk write size (trade-off between RAM and disk IO)
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
output_filename = 'data/synthetic/classification_dataset_stft_awgn.mat';
if exist(output_filename, 'file'), delete(output_filename); end
m = matfile(output_filename, 'Writable', true);

total_signals = numel(snr_levels_db) * numel(hopset_initial_states_g2) * num_signals_per_config;
total_examples = total_signals * simParams.numHops;

% Pre-allocate variables in the MAT-file
m.X_data = cell(total_examples, 1);
m.Y_data = zeros(total_examples, 1, 'uint8');
m.numHopsPerSignal = simParams.numHops;
m.stft_params = stft_params;

signal_counter = 0;
example_write_counter = 0;
chunk_idx = 0;
X_chunk = cell(chunkSize, 1);
Y_chunk = zeros(chunkSize, 1, 'uint8');

fprintf('Starting dataset generation...\n');
fprintf('Total signals: %d | Total hops: %d\n', total_signals, total_examples);

for i = 1:numel(snr_levels_db)
    snr_db = snr_levels_db(i);
    fprintf('SNR = %d dB\n', snr_db);

    for j = 1:numel(hopset_initial_states_g2)
        simParams.pnInitial2 = hopset_initial_states_g2{j};

        for k = 1:num_signals_per_config
            % --- Generate hop sequence ---
            goldCodeSequence = generateGoldCodeSequence(simParams);
            reshapedBits = reshape(goldCodeSequence, simParams.k, []);
            transposedBits = reshapedBits';
            hopIndices_col = bi2de(transposedBits, 'left-msb');
            hopIndices = hopIndices_col'; % row vector of hop indices
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

                X_chunk{chunk_idx} = single(20*log10(abs(S) + eps)); % log-magnitude
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
            if mod(signal_counter, 20) == 0
                fprintf('  Generated %d/%d signals | written %d/%d hops\n', ...
                    signal_counter, total_signals, example_write_counter, total_examples);
            end
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

fprintf('Dataset generation complete. Saved %d examples to %s\n', example_write_counter, output_filename);
