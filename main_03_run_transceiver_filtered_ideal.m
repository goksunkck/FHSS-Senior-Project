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

% --- System-Wide Parameters ---
simParams.fs = 10e6; % 10 MHz Sample Rate. This is our "world clock".

% --- Data Source & Modulation Parameters ---
simParams.numBits = 2048; % Increased for a longer simulation
simParams.M = 2; % BFSK

% --- "REALISTIC" SIGNAL PARAMETERS ---
% 1. Symbol Rate: 0.1 Msym/s (0.1 MHz). 
simParams.symbolRate = 1e5; 

% 2. Freq Separation: 1 MHz. We keep h=1 (freqSeparation = symbolRate).
simParams.freqSeparation = simParams.symbolRate; 

% 3. Samples per Symbol: This is now calculated, not defined.

simParams.samplesPerSymbol = round(simParams.fs / simParams.symbolRate);
% --- END OF REALISTIC PARAMETERS ---

% --- Hopping & PN Generator Parameters ---
simParams.k = 3; % 2^k channels
simParams.pnPoly = [5 2 0];
simParams.pnInitial = [0 0 0 0 1];
simParams.numHops = 256;

% --- Hopset Definition (THE CRITICAL FIX) ---
numChannels = 2^simParams.k; % 8 channels

% **FIX:** The spacing MUST be > signal bandwidth to prevent interference.
% Our signal BW is approx. Rs + dF = 1MHz + 1MHz = 2 MHz.
% Let's use 2 MHz as our channel spacing for a clean, non-overlapping hopset.
spacing = 2 * simParams.freqSeparation; % 2 MHz spacing

% **FIX:** Start the hopset at 0 Hz as requested.
baseFreq = 2e6; 
simParams.hopset = (0:numChannels-1) * spacing + baseFreq;

% New Hopset: [0, 2, 4, 6, 8, 10, 12, 14] MHz
% The highest signal (at 14MHz) will span ~13-15 MHz.
% This fits easily within the 30 MHz Nyquist limit (fs/2).
% --- END OF HOPSET FIX ---

% --- Intermediate Frequency (IF) Filter configuration ---
simParams.applyIFFilter = true;
% **FIX:** Set a proper order for the FIR filter in createReceiver.m
simParams.ifFilterOrder = 70;
% **FIX:** DELETE the simParams.ifFilterCutoffHz line.
% Let the receiver calculate its own optimal cutoff:
% It will be: 1.1 * (f_sep/2 + R_s) = 1.1 * (0.5M + 1M) = 1.65 MHz
% This is a perfect cutoff for our 1 MHz signal.

% --- Derived Parameters (Calculated for convenience) ---
simParams.bitsPerSymbol = log2(simParams.M);
if mod(simParams.numBits, simParams.bitsPerSymbol) ~= 0
    error('numBits must be a multiple of bitsPerSymbol (log2(M))');
end
simParams.numSymbols = simParams.numBits / simParams.bitsPerSymbol;
simParams.bitsPerFrame = simParams.numHops * simParams.k;
simParams.numModulatedSamples = simParams.numSymbols * simParams.samplesPerSymbol; 
simParams.samplesPerHop = floor(simParams.numModulatedSamples / simParams.numHops);

fprintf('Simulation parameters loaded.\n');
fprintf('Sample Rate: %.2f MHz\n', simParams.fs/1e6);
fprintf('Symbol Rate: %.2f MHz\n', simParams.symbolRate/1e6);
fprintf('Samples/Symbol: %d\n', simParams.samplesPerSymbol);
fprintf('Channel Spacing: %.2f MHz\n', spacing/1e6);
fprintf('Hopset: %s MHz\n', mat2str(simParams.hopset/1e6));


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
[S, F, T] = spectrogram(fhssSignal, hamming(windowLength_TX), overlapLength_TX, nfft_TX, simParams.fs, "centered");
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
[receivedBits, filteredDehoppedSignal, rawDehoppedSignal] = createReceiver(receivedSignal, simParams, hopFrequencies);

ax2 = subplot(2, 1, 2); % Get handle to the second axis

% ** Use the new RX parameters for high frequency resolution **
[S_dehop, F_dehop, T_dehop] = spectrogram(filteredDehoppedSignal, hamming(windowLength_RX), overlapLength_RX, nfft_RX, simParams.fs, "centered");

imagesc(ax2, T_dehop, F_dehop, 20*log10(abs(S_dehop)));
axis(ax2, 'xy');
colorbar(ax2);
title(ax2, 'Receiver''s Filtered De-hopped Signal (Spectrogram)');
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


plotSamples = 2000:6000; % Plot long enough to see bits flip
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
plot(timeVector, real(filteredDehoppedSignal(plotSamples)));
title('Receiver''s De-hopped Signal (real part)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;