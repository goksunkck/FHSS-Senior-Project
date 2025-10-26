% TEST_AWGN_CHANNEL - A test script for the FHSS/AWGN simulation
%
% This script verifies that the project's path is set up correctly and
% that the `addAWGN.m` utility function can be called.
%
% This script should be saved in the project's root directory:
%   FHSS_Pattern_Analysis/test_awgn_channel.m
%
% It replicates the logic from Chapters 1, 2, and 4.

clear; close all; clc;

% --- 0. Setup Project Path ---
% Add all folders in 'src' (e.g., src/utils, src/components)
% to the MATLAB path so we can call functions from them.
addpath(genpath('src'));
fprintf('Added src/ and all subfolders to the MATLAB path.\n');


% --- 1. Parameters (Chapters 1 & 2) ---
fs = 1e6;             % Sample rate (1 MHz)
numBits = 1024;       % Number of bits to transmit
M = 2;                % M-ary number for FSK (Binary FSK)
freqSeparation = 200e3; % Frequency separation (200 kHz)
samplesPerSymbol = 16;  % Samples per symbol
symbolRate = fs / samplesPerSymbol;

% --- 2. Data Source (Chapter 1.3) ---
messageBits = randi([0 1], numBits, 1);

% --- 3. FSK Modulation (Chapter 1.3) ---
fskModulator = comm.FSKModulator( ...
    'ModulationOrder', M, ...
    'FrequencySeparation', freqSeparation, ...
    'SamplesPerSymbol', samplesPerSymbol, ...
    'SymbolRate', symbolRate);
modulatedData = fskModulator(messageBits);

% --- 4. PN Sequence & Hopping (Chapter 2.2) ---
k = 3;                % Group bits into 3-bit words (2^3 = 8 frequencies)
hopset = (1:2^k) * 100e3 + 10e6; % e.g., 10.1MHz to 10.8MHz
numHops = 128;        % Number of hops
bitsPerFrame = numHops * k;

pnSequenceGenerator = comm.PNSequence( ...
    'Polynomial', [5 2 0], ...        % z^5 + z^2 + 1
    'InitialConditions', [0 0 0 0 1], ...
    'SamplesPerFrame', bitsPerFrame);

binarySequence = pnSequenceGenerator();
hopIndices = bi2de(reshape(binarySequence, k, []).', 'left-msb');

% --- 5. Generate Clean FHSS Signal (Chapter 1.3) ---
fhssSignal = zeros(size(modulatedData));
samplesPerHop = floor(length(modulatedData) / numHops);
t_hop = (0:samplesPerHop-1).' / fs; % Time vector for one hop

for i = 1:numHops
    % Get data segment for this hop
    idxStart = (i-1)*samplesPerHop + 1;
    idxEnd = i*samplesPerHop;
    % Handle potential mismatch if length(modulatedData) is not
    % perfectly divisible by numHops
    if i == numHops
        idxEnd = length(modulatedData);
        dataSegment = modulatedData(idxStart:idxEnd);
        t_hop_last = (0:length(dataSegment)-1).' / fs;
        hoppingCarrier = exp(1j*2*pi*hopset(hopIndices(i) + 1)*t_hop_last);
    else
        dataSegment = modulatedData(idxStart:idxEnd);
        hoppingCarrier = exp(1j*2*pi*hopset(hopIndices(i) + 1)*t_hop);
    end
    
    % Mix and store
    fhssSignal(idxStart:idxEnd) = dataSegment .* hoppingCarrier;
end

% --- 6. Apply AWGN Channel (Chapter 4.2) ---
% This is the key change: We are now calling our utility function.
snr_dB = 20; % Target SNR

% Check if the function exists (verifies our path is correct)
if exist('addAWGN', 'file')
    % Call the utility function from src/utils/addAWGN.m
    noisySignal = addAWGN(fhssSignal, snr_dB);
    fprintf('Successfully called addAWGN.m from src/utils/.\n');
    fprintf('Added AWGN to achieve a measured SNR of %d dB.\n', snr_dB);
else
    error('Could not find addAWGN.m. Make sure you are running this from the root (FHSS_Pattern_Analysis/) and src/utils/addAWGN.m exists.');
end


% --- 7. Visualization (Chapter 5) ---
windowLength = 256;
overlapLength = 192; % 75% overlap
nfft = 1024;

figure;
subplot(2, 1, 1);
spectrogram(fhssSignal, hamming(windowLength), overlapLength, nfft, fs, 'yaxis');
title('Clean FHSS Signal (Ideal Channel)');
colorbar;

subplot(2, 1, 2);
spectrogram(noisySignal, hamming(windowLength), overlapLength, nfft, fs, 'yaxis');
title(sprintf('Noisy FHSS Signal (SNR = %d dB) - Generated via addAWGN.m', snr_dB));
colorbar;
