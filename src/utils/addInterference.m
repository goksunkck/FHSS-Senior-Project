function interferedSignal = addInterference(signal, fs, hopset, jsr_dB, numJammers)
%ADDINTERFERENCE Adds narrowband interference (sinusoidal jammers) to a signal.
%   interferedSignal = addInterference(signal, fs, hopset, jsr_dB, numJammers)
%   adds 'numJammers' sinusoidal tones to the input 'signal'.
%
%   This function implements the "Narrowband Interference" scenario specified
%   in Chapter 4.2 of the research plan.
%
%   Inputs:
%       signal        - The complex input signal (e.g., fhssSignal).
%       fs            - The sample rate (in Hz).
%       hopset        - An array of the hop frequencies (in Hz). This is
%                       used to determine the frequency band of the jammers.
%       jsr_dB        - The desired Jammer-to-Signal Ratio in dB. This
%                       defines the power of each jammer relative to the
%                       signal power.
%       numJammers    - The number of jammers to add.
%
%   Outputs:
%       interferedSignal - The signal with added narrowband interference.
%
%   Ref: Chapter 4.2
%
%   See also: bandpower.

    % Ensure Signal Processing Toolbox is available for bandpower
    if ~license('test', 'Signal_Toolbox')
        error('addInterference requires the Signal Processing Toolbox for bandpower.');
    end

    % 1. Measure the power of the input signal
    % Using bandpower is more robust than simple sum(abs(signal).^2)/N
    signalPower_watts = bandpower(signal, fs, [0 fs/2]);
    
    % 2. Calculate the power for a single jammer
    jammerPower_watts = signalPower_watts * 10^(jsr_dB / 10);
    
    % 3. Calculate the required amplitude for a complex exponential jammer
    % Power of A*exp(j*w*t) is A^2
    jammerAmplitude = sqrt(jammerPower_watts);
    
    % 4. Create the time vector
    t = (0:length(signal)-1).' / fs;
    
    % 5. Initialize the output signal
    interferedSignal = signal;
    
    % 6. Define the band for placing jammers
    minFreq = min(hopset);
    maxFreq = max(hopset);
    bandWidth = maxFreq - minFreq;

    % 7. Create and add each jammer
    for i = 1:numJammers
        % Select a random frequency within the hopping band
        jammerFreq = minFreq + rand() * bandWidth;
        
        % Generate the complex sinusoidal jammer
        jammerSignal = jammerAmplitude * exp(1j*2*pi*jammerFreq*t);
        
        % Add the jammer to the signal
        interferedSignal = interferedSignal + jammerSignal;
    end

end