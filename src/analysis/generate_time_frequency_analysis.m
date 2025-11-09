function generate_time_frequency_analysis()
%GENERATE_TIME_FREQUENCY_ANALYSIS Compare WVD and SPWVD for the FHSS signal.
%   GENERATE_TIME_FREQUENCY_ANALYSIS() synthesises an FHSS waveform using the
%   project transmitter, takes a slice of it, and then
%   computes both the Wigner-Ville distribution (WVD) and the smoothed pseudo
%   Wigner-Ville distribution (SPWVD). The resulting time-frequency
%   distributions are displayed in a single figure for easy comparison.
%
%   This utility is intended for exploratory analysis and lives in the
%   analysis toolbox so that future experimentation can build on the
%   reusable WVD/SPWVD helpers provided in this directory.
%
%   See also COMPUTEWIGNERVILLEDISTRIBUTION, COMPUTESMOOTHEDPSEUDOWIGNERVILLE.

    % --- Solves the core problem of old/stale variables ---
    clear all;
    close all;
    clc;
    
    addpath('src/components');
    
    %% Simulation parameters (mirroring the main end-to-end test)
    simParams = struct();
    simParams.fs = 10e6;               % Sample rate (Hz)
    simParams.numBits = 2048;          % Bit payload length
    simParams.M = 2;                   % BFSK alphabet size
    simParams.symbolRate = 1e5;        % Symbol rate (Hz)
    simParams.freqSeparation = simParams.symbolRate; % Tone separation (Hz)
    simParams.samplesPerSymbol = round(simParams.fs / simParams.symbolRate);
    simParams.k = 3;                   % PN generator order (2^k channels)
    simParams.pnPoly = [5 2 0];
    simParams.pnInitial = [0 0 0 0 1];
    simParams.numHops = 256;
    spacing = 2 * simParams.freqSeparation;
    baseFreq = 0; % --- You changed this to 0 ---
    simParams.hopset = (0:(2^simParams.k)-1) * spacing + baseFreq;
    simParams.bitsPerSymbol = log2(simParams.M);
    simParams.numSymbols = simParams.numBits / simParams.bitsPerSymbol;
    simParams.samplesPerHop = floor((simParams.numSymbols * simParams.samplesPerSymbol) / simParams.numHops);
    simParams.numModulatedSamples = simParams.numSymbols * simParams.samplesPerSymbol;
    simParams.bitsPerFrame = simParams.numHops * simParams.k;
    
    %% Generate the FHSS signal using the project transmitter
    [fhssSignal, ~, hopFrequencies] = createTransmitter(simParams); %#ok<ASGLU>
    fprintf('Generated FHSS signal with %d samples.\n', length(fhssSignal));
    
    % --- Print the "Ground Truth" to the console ---
    fprintf('------------------------------------\n');
    fprintf('Possible Hop Frequencies (simParams.hopset):\n');
    disp(simParams.hopset' / 1e6); % Transpose and convert to MHz for easy reading
    
    fprintf('------------------------------------\n');
    fprintf('Actual Hop Sequence (first 10 hops) in MHz:\n');
    disp(hopFrequencies(1:10)' / 1e6); % Transpose and convert to MHz
    
    %% Grab a slice of the signal for manageable WVD computation
    maxAnalysisSamples = 4096; % Increased to see more hops
    if length(fhssSignal) > maxAnalysisSamples
        fprintf('Using the first %d samples for time-frequency analysis.\n', maxAnalysisSamples);
        fhssSignalAnalysis = fhssSignal(1:maxAnalysisSamples);
        fhssSignalAnalysis = fhssSignalAnalysis(:);
        analysisSampleRate = simParams.fs; % --- This is the fix for aliasing ---
    else
        fhssSignalAnalysis = fhssSignal;
        analysisSampleRate = simParams.fs;
    end
    
    %% Compute the Wigner-Ville distribution
    [timeVector, frequencyVector, wvdMatrix] = computeWignerVilleDistribution(fhssSignalAnalysis, analysisSampleRate);
    
    %% Compute the smoothed pseudo Wigner-Ville distribution (reuse the WVD)
    [~, ~, spwvdMatrix] = computeSmoothedPseudoWignerVille([], analysisSampleRate, ...
        'WVDMatrix', wvdMatrix, 'TimeVector', timeVector, 'FrequencyVector', frequencyVector);
        
    %% Convert to decibel magnitude for visualization stability
    magnitudeWVD = 20 * log10(wvdMatrix + eps);
    magnitudeSPWVD = 20 * log10(spwvdMatrix + eps);
    
    %% Plot the results
    figure('Name', 'FHSS Time-Frequency Analysis');
    
    subplot(2, 1, 1);
    imagesc(timeVector, frequencyVector / 1e6, magnitudeWVD);
    axis xy;
    colorbar;
    title('Wigner-Ville Distribution of FHSS Signal');
    xlabel('Time (s)');
    ylabel('Frequency (MHz)');
    
    subplot(2, 1, 2);
    imagesc(timeVector, frequencyVector / 1e6, magnitudeSPWVD);
    axis xy;
    colorbar;
    title('Smoothed Pseudo Wigner-Ville Distribution of FHSS Signal');
    xlabel('Time (s)');
    ylabel('Frequency (MHz)');
    
    sgtitle('FHSS Signal Time-Frequency Comparison');

end