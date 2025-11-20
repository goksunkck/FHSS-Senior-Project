% --- create_lstm_dataset_stft.m ---
% 
% Generates a dataset of FHSS signals with AWGN for use in an LSTM model.
% This version uses STFT for time-frequency representation and Gold codes
% for the hopping sequences, as per user request.
% 
% The final dataset is saved to 'data/synthetic/lstm_dataset_stft_awgn.mat'
% and contains two variables:
%   - X_data: A cell array where each cell is a [Frequency x Time] matrix
%             representing the log-magnitude of the STFT (spectrogram).
%   - Y_data: A cell array where each cell is a vector of hop indices
%             corresponding to the signal in X_data.

clear; clc; close all;

% --- 0. Setup Project Path ---
addpath(genpath('src'));
fprintf('Added src/ and all subfolders to the MATLAB path.\n');

% --- 1. Configuration ---

% SNR levels in dB
snr_levels_db = -10:2:10;

% Number of signals to generate per hopset and SNR level
num_signals_per_config = 9; % Increased to generate ~1000 signals

% Gold Code "preferred pair" polynomials of degree 5
pnPoly1 = [5 2 0];
pnPoly2 = [5 4 3 2 0];

% Define different hopsets by using different initial states for the 2nd generator.
% The register length is 5.
hopset_initial_states_g2 = { ...
    [0 0 0 0 1], ...
    [0 0 0 1 0], ...
    [0 0 1 1 1], ...
    [0 1 0 0 1], ...
    [0 1 1 0 0], ...
    [1 0 0 0 1], ...
    [1 0 1 1 0], ...
    [1 1 0 0 1], ...
    [1 1 1 0 1], ...
    [1 1 1 1 1] ...
};
num_hopsets = length(hopset_initial_states_g2);

% STFT Parameters (from main_04_run_transceiver_filtered_awgn.m)
stft_params = struct();
stft_params.windowLength = 256;
stft_params.overlapLength = 192; % 75% overlap
stft_params.nfft = 1024;


% --- 2. Dataset Initialization ---
X_data = {}; % Cell array for input sequences (spectrograms)
Y_data = {}; % Cell array for output labels (hop indices)
total_signals = length(snr_levels_db) * num_hopsets * num_signals_per_config;
signal_counter = 0;

% --- 3. Define Base Simulation Parameters (from main_04) ---
simParams = struct();
simParams.fs = 10e6;
simParams.numBits = 512; % Reduced for memory saving
simParams.M = 2;
simParams.symbolRate = simParams.fs * 0.01;
simParams.freqSeparation = simParams.symbolRate;
simParams.samplesPerSymbol = round(simParams.fs / simParams.symbolRate);
simParams.k = 3; % 2^k channels
simParams.numHops = 256;
numChannels = 2^simParams.k;
spacing = 2*simParams.freqSeparation;
baseFreq = 2e6;
simParams.hopset = (0:numChannels-1) * spacing + baseFreq;
simParams.bitsPerSymbol = log2(simParams.M);
simParams.numSymbols = simParams.numBits / simParams.bitsPerSymbol;
simParams.bitsPerFrame = simParams.numHops * simParams.k;
simParams.numModulatedSamples = simParams.numSymbols * simParams.samplesPerSymbol;
simParams.samplesPerHop = floor(simParams.numModulatedSamples / simParams.numHops);

% Gold Code specific params
simParams.pnPoly1 = pnPoly1;
simParams.pnInitial1 = [0 0 0 0 1]; % Fixed initial state for the 1st generator
simParams.pnPoly2 = pnPoly2; % Add pnPoly2 to simParams

% FSK Modulator (defined once for efficiency)
modulator = comm.FSKModulator(simParams.M, simParams.freqSeparation, ...
    'SamplesPerSymbol', simParams.samplesPerSymbol, ...
    'SymbolRate', simParams.symbolRate);

% --- 4. Data Generation Loop ---
fprintf('Starting dataset generation with STFT and Gold Codes...\n');
fprintf('Total signals to generate: %d\n', total_signals);

for i = 1:length(snr_levels_db)
    snr_db = snr_levels_db(i);
    fprintf('Processing SNR = %d dB\n', snr_db);

    for j = 1:num_hopsets
        fprintf('  Using hopset %d/%d...\n', j, num_hopsets);
        
        % Update simParams with the current Gold Code configuration
        simParams.pnInitial2 = hopset_initial_states_g2{j}; % Set pnInitial2 here
        
        for k = 1:num_signals_per_config
            % --- A. Generate Hop Sequence using Gold Code ---
            % 1. Generate the Gold code bit sequence
            goldCodeSequence = generateGoldCodeSequence(simParams);
            
            % 2. Replicate logic from generateHopIndices.m
            % Reshape the bit stream into a matrix where each column has k bits
            reshapedBits = reshape(goldCodeSequence, simParams.k, []);
            % Transpose so each row is a k-bit word, and convert binary to decimal
            hopIndices = bi2de(reshapedBits', 'left-msb')'; % hopIndices is now a row vector

            % Map indices to frequencies for signal generation
            hopFrequencies = simParams.hopset(hopIndices + 1); % +1 for 1-based MATLAB indexing

            % --- B. Generate FHSS signal (replicating createTransmitter logic) ---
            % 1. Generate message bits and modulate
            messageBits = randi([0 1], simParams.numBits, 1);
            modulatedData = modulator(messageBits);
            
            % 2. Mix with hopping carrier
            fhssSignal = zeros(simParams.numModulatedSamples, 1);
            for hop_idx = 1:simParams.numHops
                startIndex = (hop_idx-1)*simParams.samplesPerHop + 1;
                stopIndex = hop_idx*simParams.samplesPerHop;
                dataSegment = modulatedData(startIndex:stopIndex);
                
                t = (0:simParams.samplesPerHop-1)' / simParams.fs;
                hoppingCarrier = exp(1j*2*pi * hopFrequencies(hop_idx) * t);
                
                fhssSignal(startIndex:stopIndex) = dataSegment .* hoppingCarrier;
            end
            
            % --- C. Apply AWGN Channel ---
            receivedSignal = addAWGN(fhssSignal, snr_db);
            
            % --- D. Perform STFT Analysis ---
            [S, ~, ~] = spectrogram(receivedSignal, ...
                                    hamming(stft_params.windowLength), ...
                                    stft_params.overlapLength, ...
                                    stft_params.nfft, ...
                                    simParams.fs, ...
                                    "centered");
            
            % Use log-magnitude of the STFT as the feature matrix
            stft_log_magnitude = 20*log10(abs(S));

            % --- E. Store Data ---
            % X_data will be the spectrogram matrix
            X_data{end+1} = stft_log_magnitude;
            % Y_data will be the vector of hop indices
            Y_data{end+1} = hopIndices;
            
            signal_counter = signal_counter + 1;
            if mod(signal_counter, 50) == 0
                fprintf('    Generated %d/%d signals...\n', signal_counter, total_signals);
            end
        end
    end
end

% --- 5. Save Dataset ---
output_filename = 'data/synthetic/lstm_dataset_stft_awgn.mat';
fprintf('Saving dataset to %s...\n', output_filename);
save(output_filename, 'X_data', 'Y_data', '-v7.3');

fprintf('Dataset generation complete. Saved %d signals.\n', signal_counter);
