function hopIndices = generateHopIndices(pnGen, k)
    % Generates a frame of hop indices from a PN generator object.
    % Based on Chapter 2.2, Code Example
    %
    % pnGen: A comm.PNSequence System Object
    % k:     Number of bits per hop (e.g., k=3 for 2^3=8 frequencies)
    
    % Generate the binary sequence
    binarySequence = pnGen();
    
    % Reshape the binary sequence into k-bit words and convert to decimal indices
    numWords = length(binarySequence) / k;
    hopIndices = bi2de(reshape(binarySequence, k, numWords)', 'left-msb'); 
end