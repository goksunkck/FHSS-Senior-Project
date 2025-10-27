function [receivedBits, filteredDehoppedSignal, rawDehoppedSignal] = createReceiver(receivedSignal, simParams, hopFrequencies)
    % Implements the FHSS Receiver from Chapter 1.4
    
    %% 1. De-hopping (Ch 1.4, Step 1)
    %
    dehoppedSignal = zeros(length(receivedSignal), 1);
    
    for i = 1:simParams.numHops
        % Get the segment of the received signal for this hop
        startIndex = (i-1)*simParams.samplesPerHop + 1;
        stopIndex = i*simParams.samplesPerHop;
        signalSegment = receivedSignal(startIndex:stopIndex);
        
        % Generate the local hopping carrier (must be identical to TX)
        t = (0:simParams.samplesPerHop-1)' / simParams.fs;
        localHoppingCarrier = exp(1j*2*pi * hopFrequencies(i) * t);
        
        % De-hop by multiplying with the conjugate
        dehoppedSignal(startIndex:stopIndex) = signalSegment .* conj(localHoppingCarrier);
    end
    rawDehoppedSignal = dehoppedSignal; % Save the raw dehopped signal for output
    
    %% 2. IF Channel Filtering (Ch 1.4, Step 2)
    %
    % Apply an optional channel filter (low-pass) prior to the
    % demodulator. This is crucial for removing out-of-band noise and 
    % interference that was not removed by the RF front-end.
    
    filteredDehoppedSignal = rawDehoppedSignal; % Default to raw if no filter
    applyFilter = isfield(simParams, 'applyIFFilter') && simParams.applyIFFilter;
    
if applyFilter
        filterType = 'lowpass';
        if isfield(simParams, 'ifFilterType')
            filterType = lower(simParams.ifFilterType);
        end
        
        switch filterType
            case 'lowpass'
                % --- Linear Phase FIR (Correct for Communications) ---
                % This filter provides a sharp cutoff (good spectrogram)
                % without introducing phase distortion (good BER).

                % Set a default FIR order (higher is sharper)
                if isfield(simParams, 'ifFilterOrder')
                    filterOrder = simParams.ifFilterOrder;
                else
                    filterOrder = 70; % Default for FIR
                end
                
                % Use a dynamic cutoff if one isn't provided
                if isfield(simParams, 'ifFilterCutoffHz')
                    cutoffHz = simParams.ifFilterCutoffHz;
                else
                    % Optimal cutoff: Pass the tone + main lobe
                    maxSignalFreq = (simParams.freqSeparation / 2) + simParams.symbolRate;
                    guardBandFactor = 1.1; 
                    cutoffHz = maxSignalFreq * guardBandFactor;
                end

                % Cache the filter object for speed
                persistent cachedFilterObject cachedFs cachedCutoff cachedOrder
                if isempty(cachedFilterObject) || cachedFs ~= simParams.fs || ...
                        cachedCutoff ~= cutoffHz || cachedOrder ~= filterOrder
                    
                    % Use 'designfilt' to create a linear-phase FIR filter
                    % 'kaiserwin' gives a good tradeoff between sharpness
                    % and stopband attenuation.
                    cachedFilterObject = designfilt('lowpassfir', ...
                        'FilterOrder', filterOrder, ...
                        'CutoffFrequency', cutoffHz, ...
                        'SampleRate', simParams.fs, ...
                        'Window', 'kaiser');
                    
                    cachedFs = simParams.fs;
                    cachedCutoff = cutoffHz;
                    cachedOrder = filterOrder;
                end
                
                % Use 'filtfilt' for zero-phase filtering.
                % This is ideal for simulation as it introduces ZERO
                % phase distortion, solving the ISI/BER problem.
                filteredDehoppedSignal = filtfilt(cachedFilterObject, rawDehoppedSignal);

            otherwise
                warning('Unknown IF filter type "%s". Skipping IF filtering.', filterType);
        end
    end
    
    % NOTE: The critical bug 'filteredDehoppedSignal = dehoppedSignal;'
    % has been removed from here.


    %% 3. Demodulation (Ch 1.4, Step 2)
    %
    % The FSKDemodulator object contains the optimal internal matched
    % filters (or filter bank) needed to demodulate the BFSK signal.
    demodulator = comm.FSKDemodulator(simParams.M, simParams.freqSeparation, ...
        'SamplesPerSymbol', simParams.samplesPerSymbol, ...
        'SymbolRate', simParams.symbolRate); 
    
    % Pass the *channel-filtered* dehopped signal to the demodulator
    receivedBits = demodulator(filteredDehoppedSignal); 
end