% -----------------------------------------------------------------
% main_01_run_transceiver.m
%
% This script runs a complete, end-to-end FHSS transceiver simulation
% under ideal (noise-free) conditions. It verifies that the
% transmitter and receiver logic is working correctly before
% moving on to dataset generation.
%
% Implements logic from Chapters 1 & 2.
% -----------------------------------------------------------------
clear; clc; close all; clear functions;
addpath('src/components');
%% 1. Define Simulation Parameters
%
simParams = struct();

% --- FIXES ARE HERE ---
% System-Wide Parameters
simParams.fs = 60e6; % Sample Rate (Hz) <-- INCREASED from 1e6 to 60e6

% Data Source & Modulation Parameters (Chapter 1)
simParams.numBits = 1024;
simParams.M = 2;
simParams.freqSeparation = 250e3;

% ** We must increase SamplesPerSymbol to keep the SymbolRate the same **
old_fs = 1e6;
old_sps = 16;
old_symbolRate = old_fs / old_sps; % This was 62,500 symbols/sec

% New SamplesPerSymbol to maintain the same SymbolRate
simParams.samplesPerSymbol = round(simParams.fs / old_symbolRate); % <-- ADJUSTED
% --- END OF FIXES ---

% Hopping & PN Generator Parameters (Chapters 1 & 2)
simParams.k = 3;
simParams.pnPoly = [5 2 0];
simParams.pnInitial = [0 0 0 0 1];
simParams.numHops = 256;
simParams.hopset = [10e6, 12e6, 14e6, 16e6, 18e6, 20e6, 22e6, 24e6];

% Derived Parameters (Calculated for convenience)
simParams.bitsPerSymbol = log2(simParams.M); % log2(M) bits/symbol
if mod(simParams.numBits, simParams.bitsPerSymbol) ~= 0
    error('numBits must be a multiple of bitsPerSymbol (log2(M))');
end
simParams.numSymbols = simParams.numBits / simParams.bitsPerSymbol;
simParams.symbolRate = simParams.fs / simParams.samplesPerSymbol;
simParams.bitsPerFrame = simParams.numHops * simParams.k;

% **CRITICAL FIX: Base this on numSymbols, not numBits**
simParams.numModulatedSamples = simParams.numSymbols * simParams.samplesPerSymbol; 
simParams.samplesPerHop = floor(simParams.numModulatedSamples / simParams.numHops);

fprintf('Simulation parameters loaded.\n');
fprintf('New Sample Rate: %.2f MHz\n', simParams.fs/1e6);
fprintf('New Symbol Rate: %.2f kHz\n', simParams.symbolRate/1e3);
fprintf('New Samples/Symbol: %d\n', simParams.samplesPerSymbol);

%% 2. Run Transmitter
%
% Generate the message bits and the final FHSS signal.
% We also get the 'hopFrequencies' sequence for the receiver.
[fhssSignal, messageBits, hopFrequencies, modulatedData] = createTransmitter(simParams);

fprintf('Transmitter complete. Generated %d samples.\n', length(fhssSignal));

%% 3. Ideal Channel
%
% For this initial test, we assume a perfect, noise-free channel.
% The received signal is identical to the transmitted signal.
receivedSignal = fhssSignal;

%% 4. Run Receiver
%
% The receiver uses the same parameters to de-hop and demodulate.
% We pass the 'hopFrequencies' to simulate perfect synchronization.
receivedBits = createReceiver(receivedSignal, simParams, hopFrequencies);

fprintf('Receiver complete. Recovered %d bits.\n', length(receivedBits));

%% 5. Performance Validation
%
% Compare the original transmitted bits to the recovered bits
% to calculate the Bit Error Rate (BER).
[numErrors, ber] = biterr(messageBits, receivedBits); 

fprintf('\n--- Simulation Results ---\n');
fprintf('Total Bits Transmitted: %d\n', simParams.numBits);
fprintf('Total Bit Errors: %d\n', numErrors);
fprintf('Bit Error Rate (BER): %f\n', ber);

if ber == 0
    fprintf('\nSuccess! The transceiver link is working perfectly.\n');
else
    fprintf('\nWarning: Errors detected. Check transceiver logic.\n');
end
%% 6. Plot Signals
%
% Let's visualize the signals using the STFT (Spectrogram)
%
% ** FIX: Increase 'windowLength' for better frequency resolution.

% --- STFT Parameters for TX Plot (Good time resolution) ---
windowLength_TX = 256;
overlapLength_TX = 192; % 75% overlap
nfft_TX = 1024; % Zero-padding

% --- STFT Parameters for RX Plot (Good frequency resolution) ---
windowLength_RX = 4096; % <-- INCREASED for sharp frequencies
overlapLength_RX = 3072; % <-- INCREASED (75% overlap)
nfft_RX = 4096; % <-- INCREASED (must be >= windowLength)

% --- 1. Plot the Transmitted FHSS Signal ---
figure;
ax1 = subplot(2, 1, 1); % Get handle to the first axis
[S, F, T] = spectrogram(fhssSignal, hamming(windowLength_TX), overlapLength_TX, nfft_TX, simParams.fs);
imagesc(ax1, T, F, 20*log10(abs(S)));
axis(ax1, 'xy');
colorbar(ax1);
title(ax1, 'Transmitted FHSS Signal (Spectrogram)');
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Frequency (Hz)');
ylim_min = min(simParams.hopset) - 2*simParams.freqSeparation;
ylim_max = max(simParams.hopset) + 2*simParams.freqSeparation;
ylim(ax1, [ylim_min, ylim_max]);

% --- 2. Plot the De-hopped Signal at the Receiver ---
[receivedBits, dehoppedSignal] = createReceiver(receivedSignal, simParams, hopFrequencies);

ax2 = subplot(2, 1, 2); % Get handle to the second axis

% ** Use the new RX parameters for high frequency resolution **
[S_dehop, F_dehop, T_dehop] = spectrogram(dehoppedSignal, hamming(windowLength_RX), overlapLength_RX, nfft_RX, simParams.fs, "centered");

imagesc(ax2, T_dehop, F_dehop, 20*log10(abs(S_dehop)));
axis(ax2, 'xy');
colorbar(ax2);
title(ax2, 'Receiver''s De-hopped Signal (Spectrogram)');
xlabel(ax2, 'Time (s)');
ylabel(ax2, 'Frequency (Hz)');

% ** This Y-limit will now work correctly **
ylim(ax2, [-(simParams.freqSeparation), (simParams.freqSeparation)]); 


%% 7. Plot Bits
%
% Let's plot the original and recovered bits to see the BER of 0.

figure;
subplot(2, 1, 1);
stem(messageBits, 'b.');
title('Transmitted Bits');
xlabel('Bit Index');
ylabel('Value (0 or 1)');
% Only plot the first ~200 bits so we can see them clearly
xlim([0, 200]);
ylim([-0.2, 1.2]);

subplot(2, 1, 2);
stem(receivedBits, 'r.');
title('Received Bits');
xlabel('Bit Index');
ylabel('Value (0 or 1)');
% Plot the same 200 bits for comparison
xlim([0, 200]);
ylim([-0.2, 1.2]);

%% 8. Plot Waveforms
%
% Plot the real part of the signals in the time domain.


plotSamples = 20000:60000; % Plot long enough to see bits flip
timeVector = (plotSamples - 1) / simParams.fs; % Time axis in seconds

figure;

% --- 1. Baseband FSK Signal (before hopping) ---
subplot(3, 1, 1);
plot(timeVector, real(modulatedData(plotSamples)));
title('Baseband FSK Signal (real part)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- 2. Transmitted FHSS Signal (after hopping) ---
subplot(3, 1, 2);
plot(timeVector, real(fhssSignal(plotSamples)));
title('Transmitted FHSS Signal (real part)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% --- 3. De-hopped Signal (at receiver) ---
subplot(3, 1, 3);
plot(timeVector, real(dehoppedSignal(plotSamples)));
title('Receiver''s De-hopped Signal (real part)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;