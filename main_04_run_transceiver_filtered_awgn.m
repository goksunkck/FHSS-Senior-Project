% -----------------------------------------------------------------
% main_02_run_transceiver_awgn.m
%
% This script runs a complete, end-to-end FHSS transceiver simulation
% under AWGN (noisy) conditions. It verifies the transceiver's
% performance (BER) in a non-ideal channel.
%
% This is based on main_01_run_transceiver.m, but replaces the
% "Ideal Channel" section with a call to the addAWGN utility.
%
% Implements logic from Chapters 1, 2, 3, and 4.
% -----------------------------------------------------------------
clear; clc; close all; clear functions;

% --- 0. Setup Project Path ---
% Add all folders in 'src' (e.g., src/utils, src/components)
% to the MATLAB path so we can call functions from them.
addpath(genpath('src'));
fprintf('Added src/ and all subfolders to the MATLAB path.\n');

%% 1. Define Simulation Parameters
%
simParams = struct();
% System-Wide Parameters
simParams.fs = 10e6; % Sample Rate (Hz)
simParams.snr_dB = 0; % Signal-to-Noise Ratio in dB (Adjusted for new PG)
% Data Source & Modulation Parameters (Chapter 1)
simParams.numBits = 4096;
simParams.M = 2;

% --- WIDEBAND PARAMETERS (Rs = 0.1 * fs) ---
% Set symbol rate to 10% of the sample rate (Rs = 0.1*fs)
simParams.symbolRate = simParams.fs * 0.01; % 0.1 MHz

% We'll use non-coherent FSK with h=1 (robust)
% So, freqSeparation = symbolRate
simParams.freqSeparation = simParams.symbolRate; 

% Calculate SamplesPerSymbol based on the new, clean parameters.
simParams.samplesPerSymbol = round(simParams.fs / simParams.symbolRate); 
% --- END OF WIDEBAND PARAMETERS ---

% Hopping & PN Generator Parameters (Chapters 1 & 2)
simParams.k = 3; % 2^k channels 
simParams.pnPoly = [5 2 0];
simParams.pnInitial = [0 0 0 0 1];
simParams.numHops = 256;
% Programmatically create a hopset. Spacing must match freqSeparation
% to avoid adjacent channel interference.
numChannels = 2^simParams.k; % 8 channels
spacing = 2*simParams.freqSeparation; % 2*freq sep spacings
baseFreq = 2e6;
simParams.hopset = (0:numChannels-1) * spacing + baseFreq;
% --- END OF HOPSET UPDATE ---

% Intermediate Frequency (IF) Filter configuration
simParams.applyIFFilter = false;
simParams.ifFilterCutoffHz = simParams.symbolRate; % Low-pass cutoff (Hz)
simParams.ifFilterOrder = 70;        

% Derived Parameters (Calculated for convenience)
simParams.bitsPerSymbol = log2(simParams.M); % log2(M) bits/symbol
if mod(simParams.numBits, simParams.bitsPerSymbol) ~= 0
    error('numBits must be a multiple of bitsPerSymbol (log2(M))');
end
simParams.numSymbols = simParams.numBits / simParams.bitsPerSymbol;
% simParams.symbolRate = simParams.fs / simParams.samplesPerSymbol; % This is now defined above
simParams.bitsPerFrame = simParams.numHops * simParams.k;
simParams.numModulatedSamples = simParams.numSymbols * simParams.samplesPerSymbol; 
simParams.samplesPerHop = floor(simParams.numModulatedSamples / simParams.numHops);
fprintf('Simulation parameters loaded.\n');
fprintf('Sample Rate: %.2f MHz\n', simParams.fs/1e6);
fprintf('Symbol Rate: %.2f MHz\n', simParams.symbolRate/1e6); % Changed to MHz
fprintf('Samples/Symbol: %d\n', simParams.samplesPerSymbol);
fprintf('Target SNR: %d dB\n', simParams.snr_dB);
fprintf('Hopset: %s MHz\n', mat2str(simParams.hopset/1e6));
%% 2. Run Transmitter
%
% Generate the message bits and the final FHSS signal.
% We also get the 'hopFrequencies' sequence for the receiver.
[fhssSignal, messageBits, hopFrequencies, modulatedData] = createTransmitter(simParams);
fprintf('Transmitter complete. Generated %d samples.\n', length(fhssSignal));

%% 3. AWGN Channel (Chapter 4.2)
%
% This replaces the "Ideal Channel" section.
% We call our utility function from src/utils/ to add noise.
fprintf('Applying AWGN channel (SNR = %d dB)...\n', simParams.snr_dB);

% Check if the function exists (verifies our path is correct)
if ~exist('addAWGN', 'file')
    error('Could not find addAWGN.m. Make sure src/utils/ is on the path.');
end

% Call the utility function from src/utils/addAWGN.m
receivedSignal = addAWGN(fhssSignal, simParams.snr_dB);


%% 4. Run Receiver
%
% The receiver uses the same parameters to de-hop and demodulate.
% We pass the 'hopFrequencies' to simulate perfect synchronization.
% (Note: createReceiver is called again in plotting section
%  to get the intermediate dehoppedSignal)
[receivedBits, ~] = createReceiver(receivedSignal, simParams, hopFrequencies);
fprintf('Receiver complete. Recovered %d bits.\n', length(receivedBits));

%% 5. Performance Validation
%
% Compare the original transmitted bits to the recovered bits
% to calculate the Bit Error Rate (BER).
[numErrors, ber] = biterr(messageBits, receivedBits); 

fprintf('\n--- Simulation Results ---\n');
fprintf('Target SNR: \t\t%d dB\n', simParams.snr_dB);
fprintf('Total Bits Transmitted: %d\n', simParams.numBits);
fprintf('Total Bit Errors: \t%d\n', numErrors);
fprintf('Bit Error Rate (BER): \t%f\n', ber);

if ber == 0
    fprintf('\nSuccess! The transceiver link is working perfectly.\n');
else
    fprintf('\nNote: Errors detected, as expected in a noisy channel.\n');
end

%% 6. Plot Signals
%
% Let's visualize the signals using the STFT (Spectrogram)
%
% STFT Parameters
windowLength_TX = 256;
overlapLength_TX = 192; % 75% overlap
nfft_TX = 1024; 
windowLength_RX = 4096; 
overlapLength_RX = 3072; % 75% overlap
nfft_RX = 4096; 

% --- 1. Plot the Transmitted FHSS Signal (Noisy) ---
%    (We plot the receivedSignal to see the noise)
figure;
ax1 = subplot(2, 1, 1); 
[S, F, T] = spectrogram(receivedSignal, hamming(windowLength_TX), overlapLength_TX, nfft_TX, simParams.fs, "centered");
imagesc(ax1, T, F, 20*log10(abs(S)));
axis(ax1, 'xy');
colorbar(ax1);
title(ax1, sprintf('Received FHSS Signal (Spectrogram), SNR = %d dB', simParams.snr_dB));
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Frequency (Hz)');
ylim_min = min(simParams.hopset) - 2*simParams.freqSeparation;
ylim_max = max(simParams.hopset) + 2*simParams.freqSeparation;
ylim(ax1, [ylim_min, ylim_max]);

% --- 2. Plot the De-hopped Signal at the Receiver ---
%    (This signal should be noisy but centered at baseband)
[receivedBits, filteredDehoppedSignal, rawDehoppedSignal] = createReceiver(receivedSignal, simParams, hopFrequencies);
ax2 = subplot(2, 1, 2); 
[S_dehop, F_dehop, T_dehop] = spectrogram(filteredDehoppedSignal, hamming(windowLength_RX), overlapLength_RX, nfft_RX, simParams.fs, "centered");
imagesc(ax2, T_dehop, F_dehop, 20*log10(abs(S_dehop)));
axis(ax2, 'xy');
colorbar(ax2);
title(ax2, 'Receiver''s Filtered De-hopped Signal (Spectrogram)');
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Frequency (Hz)');
% Adjust Y-limits for new 6MHz freq separation (h=1)
% The dehopped signal will be at 0 and freqSeparation
ylim(ax2, [-(simParams.freqSeparation), (simParams.freqSeparation)]); 



%% 7. Plot Bits
%
% Let's plot the original and recovered bits to see the errors
figure;
subplot(2, 1, 1);
stem(messageBits, 'b.');
title('Transmitted Bits');
xlabel('Bit Index');
ylabel('Value (0 or 1)');
xlim([0, 200]);
ylim([-0.2, 1.2]);

subplot(2, 1, 2);
stem(receivedBits, 'r.');
title(sprintf('Received Bits (BER = %.4f)', ber));
xlabel('Bit Index');
ylabel('Value (0 or 1)');
xlim([0, 200]);
ylim([-0.2, 1.2]);

%% 8. Plot Waveforms
%
% Plot the real part of the signals in the time domain.
% **FIX:** The old range (20000:60000) is outside the new array bounds.
% The total signal length is numModulatedSamples = 10240.
% We will use a new, valid range to zoom in.
plotSamples = 2000:6000; 
timeVector = (plotSamples - 1) / simParams.fs; 
figure;

% --- 1. Baseband FSK Signal (before hopping) ---
subplot(3, 1, 1);
plot(timeVector, real(modulatedData(plotSamples)));
title('Baseband FSK Signal (Clean)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- 2. Received FHSS Signal (Noisy) ---
subplot(3, 1, 2);
plot(timeVector, real(receivedSignal(plotSamples)));
title(sprintf('Received FHSS Signal (Noisy, SNR = %d dB)', simParams.snr_dB));
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- 3. De-hopped Signal (at receiver) ---
subplot(3, 1, 3);
plot(timeVector, real(filteredDehoppedSignal(plotSamples)));
title('Receiver''s Filtered De-hopped Signal (Noisy)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;






