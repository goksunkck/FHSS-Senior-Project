function [fhssSignal, messageBits, hopFrequencies, modulatedData] = createTransmitter(simParams)
    % Implements the FHSS Transmitter from Chapter 1.3
    %
    % simParams: A struct containing all simulation parameters
    
    %% 1. Data Source (Ch 1.3, Step 1)
    %
    messageBits = randi([0 1], simParams.numBits, 1);

  
    %% 2. Modulation (Ch 1.3, Step 2)
    %
    modulator = comm.FSKModulator(simParams.M, simParams.freqSeparation, ...
        'SamplesPerSymbol', simParams.samplesPerSymbol, ...
        'SymbolRate', simParams.symbolRate); 
    
    modulatedData = modulator(messageBits); 

    %% 3. Frequency Synthesizer and Hopping (Ch 1.3, Step 3 & Ch 2)
    %
    % Create the PN generator
    pnGen = comm.PNSequence('Polynomial', simParams.pnPoly, ...
                           'InitialConditions', simParams.pnInitial, ...
                           'SamplesPerFrame', simParams.bitsPerFrame); 
                       
    % Generate hop indices
    hopIndices = generateHopIndices(pnGen, simParams.k);
    
    % Map indices to frequencies
    hopFrequencies = simParams.hopset(hopIndices + 1); % +1 for 1-based MATLAB indexing

    %% 4. Final Signal Generation (Ch 1.3, Step 4)
    %
    fhssSignal = zeros(simParams.numModulatedSamples, 1);
    
    for i = 1:simParams.numHops
        % Get the segment of modulated data for this hop
        startIndex = (i-1)*simParams.samplesPerHop + 1;
        stopIndex = i*simParams.samplesPerHop;
        dataSegment = modulatedData(startIndex:stopIndex); 
        
        % Generate the hopping carrier for this hop
        t = (0:simParams.samplesPerHop-1)' / simParams.fs;
        hoppingCarrier = exp(1j*2*pi * hopFrequencies(i) * t); 
        
        % Mix and store
        fhssSignal(startIndex:stopIndex) = dataSegment .* hoppingCarrier; 
    end
end