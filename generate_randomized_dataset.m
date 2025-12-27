% --- generate_randomized_dataset.m ---
%
% Generates a dataset with DEGREE 10 polynomials and FILTERED RANDOMIZED SEEDS.
%
% Output: data/synthetic/classification_dataset_stft_random_deg10.mat

clear; clc; close all;

addpath(genpath('src'));
fprintf('Added src/ and subfolders to the MATLAB path.\n');

% --- Configuration ---
snr_levels_db = -10:2:10; % 11 levels
num_signals_per_snr = 30; % Total signals = 11 * 30 = 330

% Degree 10 Polynomials (Primitives)
% Poly1: x^10 + x^3 + 1 -> [10 3 0]
% Poly2: x^10 + x^7 + 1 -> [10 7 0] (Commonly used, though strict preferred pair check might vary)
pnPoly1 = [10 3 0];
pnPoly2 = [10 7 0];

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

% Update PN Params for Degree 10
simParams.pnPoly1 = pnPoly1;
% Initial state must be length 10
simParams.pnInitial1 = [0 0 0 0 0 0 0 0 0 1];
simParams.pnPoly2 = pnPoly2;

% FSK modulator
modulator = comm.FSKModulator(simParams.M, simParams.freqSeparation, ...
    'SamplesPerSymbol', simParams.samplesPerSymbol, ...
    'SymbolRate', simParams.symbolRate);

% --- MAT-file initialization ---
output_filename = 'data/synthetic/classification_dataset_stft_random_deg10.mat';
if exist(output_filename, 'file'), delete(output_filename); end
m = matfile(output_filename, 'Writable', true);

total_signals = numel(snr_levels_db) * num_signals_per_snr;
total_examples = total_signals * simParams.numHops;

% Pre-allocate variables
m.X_data = cell(total_examples, 1);
m.Y_data = zeros(total_examples, 1, 'uint8');
m.numHopsPerSignal = simParams.numHops;
m.stft_params = stft_params;
% Store degree info
m.polynomial_degree = 10;

signal_counter = 0;
example_write_counter = 0;
chunk_idx = 0;
X_chunk = cell(chunkSize, 1);
Y_chunk = zeros(chunkSize, 1, 'uint8');

% Global Min/Max tracking
global_min = Inf;
global_max = -Inf;

fprintf('Starting RANDOMIZED dataset generation (Degree 10)...\n');
fprintf('Total signals: %d | Total hops: %d\n', total_signals, total_examples);

for i = 1:numel(snr_levels_db)
    snr_db = snr_levels_db(i);
    fprintf('SNR = %d dB\n', snr_db);

    for k = 1:num_signals_per_snr
        % --- RANDOMIZED SEED GENERATION ---
        % Generate random vector of length 10
        rand_seed = randi([0 1], 1, 10);

        % Ensure not all zeros (rule of LFSR)
        while all(rand_seed == 0)
            rand_seed = randi([0 1], 1, 10);
        end

        simParams.pnInitial2 = rand_seed;

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

            % Log-Magnitude
            val_stft = single(20*log10(abs(S) + eps));

            X_chunk{chunk_idx} = val_stft;
            Y_chunk(chunk_idx) = uint8(hopIndices(hop_idx));

            % Track min/max
            current_min = min(val_stft, [], 'all');
            current_max = max(val_stft, [], 'all');
            if current_min < global_min, global_min = current_min; end
            if current_max > global_max, global_max = current_max; end

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
        if mod(signal_counter, 50) == 0
            fprintf('  Generated %d/%d signals...\n', signal_counter, total_signals);
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

m.global_min = global_min;
m.global_max = global_max;
fprintf('\nGlobal Stats: Min = %.4f | Max = %.4f\n', global_min, global_max);

fprintf('Dataset generation complete. Saved %d examples to %s\n', example_write_counter, output_filename);
