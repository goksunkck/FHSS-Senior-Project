clear all; clc; close all;
addpath(genpath('src'));
%% Monte Carlo Simulation
% ---------- 1) common sim params ----------
baseParams = struct();
baseParams.fs               = 10e6;
baseParams.numBits          = 1024;
baseParams.M                = 2;
baseParams.symbolRate       = baseParams.fs * 0.01;
baseParams.freqSeparation   = baseParams.symbolRate;
baseParams.samplesPerSymbol = round(baseParams.fs / baseParams.symbolRate);
baseParams.k                = 3;
baseParams.pnPoly           = [5 2 0];
baseParams.pnInitial        = [0 0 0 0 1];
baseParams.numHops          = 256;

numChannels = 2^baseParams.k;
spacing    = 2 * baseParams.freqSeparation;
baseFreq   = 2e6;
baseParams.hopset = (0:numChannels-1) * spacing + baseFreq;

baseParams.applyIFFilter   = true;
baseParams.ifFilterCutoffHz= baseParams.symbolRate;
baseParams.ifFilterOrder   = 70;

baseParams.bitsPerSymbol     = log2(baseParams.M);
baseParams.numSymbols        = baseParams.numBits / baseParams.bitsPerSymbol;
baseParams.bitsPerFrame      = baseParams.numHops * baseParams.k;
baseParams.numModulatedSamples = baseParams.numSymbols * baseParams.samplesPerSymbol;
baseParams.samplesPerHop     = floor(baseParams.numModulatedSamples / baseParams.numHops);

% ---------- 2) Monte Carlo settings ----------
snrVec_dB   = -10:5:10;     % SNR points
targetErrs  = 1000;        % stop per SNR when we saw >= this many bit errors
maxFrames   = 5000;        % safety cap

berMC = zeros(size(snrVec_dB));

for iS = 1:numel(snrVec_dB)
    snrNow = snrVec_dB(iS);
    fprintf('\n--- SNR = %d dB ---\n', snrNow);

    totalErrs = 0;
    totalBits = 0;

    frameCnt  = 0;
    while (totalErrs < targetErrs) && (frameCnt < maxFrames)
        frameCnt = frameCnt + 1;

        % 2.a copy base params and set SNR for this run
        simParams = baseParams;
        simParams.snr_dB = snrNow;

        % 2.b transmitter
        [txSig, txBits, hopFreqs, ~] = createTransmitter(simParams);

        % 2.c channel
        rxSig = addAWGN(txSig, simParams.snr_dB);

        % 2.d receiver
        [rxBits, ~] = createReceiver(rxSig, simParams, hopFreqs);

        % 2.e count errors
        [nErr, ~] = biterr(txBits, rxBits);

        totalErrs = totalErrs + nErr;
        totalBits = totalBits + numel(txBits);

        fprintf('frame %3d: errs=%3d  totErr=%4d  totBits=%6d\n', ...
                frameCnt, nErr, totalErrs, totalBits);
    end

    berMC(iS) = totalErrs / totalBits;
    fprintf('BER @ %d dB = %g (from %d bits)\n', snrNow, berMC(iS), totalBits);
end

% ---------- 3) plot ----------
figure;
hold on;
set(gca, 'YScale', 'log') % semilogy and hold on works weirdly so we need to force yscale to be log like this
semilogy(snrVec_dB, berMC, 'o-');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('FHSS link BER (Monte Carlo)');


