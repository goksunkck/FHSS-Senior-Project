function [receivedBits, dehoppedSignal] = createReceiver(receivedSignal, simParams, hopFrequencies)
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

    %% 2. IF Filtering (Ch 1.4, Step 2)
    %
    % ** Per your request, no IF filter is applied **
    % The dehopped signal is passed directly to the demodulator.
    filteredDehoppedSignal = dehoppedSignal;

    %% 3. Demodulation (Ch 1.4, Step 2)
    %
    demodulator = comm.FSKDemodulator(simParams.M, simParams.freqSeparation, ...
        'SamplesPerSymbol', simParams.samplesPerSymbol, ...
        'SymbolRate', simParams.symbolRate); 
    
    % The demodulator now gets the raw, unfiltered dehopped signal
    receivedBits = demodulator(filteredDehoppedSignal); 
end