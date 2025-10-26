function noisySignal = addAWGN(signal, snr_dB)
%ADD_AWGN Adds Additive White Gaussian Noise to a signal.
%   noisySignal = addAWGN(signal, snr_dB) adds AWGN to the input 'signal'
%   to achieve the target 'snr_dB' (in decibels).
%
%   This function is a wrapper for MATLAB's 'awgn' function[cite: 204].
%   It specifically uses the 'measured' option as specified in the 
%   research plan  to ensure the SNR is calculated based 
%   on the measured power of the input signal, which is
%   critical for dataset generation.
%
%   Inputs:
%       signal   - The complex or real input signal (e.g., fhssSignal).
%       snr_dB   - The desired Signal-to-Noise Ratio in dB.
%
%   Outputs:
%       noisySignal - The signal with added AWGN.
%
%   Ref: Chapter 4.2 , Table 1 
%
%   See also: awgn.

    % Add AWGN using the "measured" flag to ensure the SNR is
    % relative to the actual signal power. 
    noisySignal = awgn(signal, snr_dB, 'measured'); 

end