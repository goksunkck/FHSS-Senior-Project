function [timeVector, frequencyVector, wvdMatrix] = computeWignerVilleDistribution(signal, sampleRate)
%COMPUTEWIGNERVILLEDISTRIBUTION Compute the Wigner-Ville distribution.
%   [T, F, W] = COMPUTEWIGNERVILLEDISTRIBUTION(X, FS) returns the time
%   vector T, frequency vector F, and the Wigner-Ville distribution matrix W
%   for the input signal X sampled at FS Hz. The implementation uses the
%   analytic representation of the signal and evaluates the discrete
%   Wigner-Ville distribution through the lag-domain FFT approach. Only the
%   non-negative frequency portion of the distribution is returned.
%
%   Inputs:
%       signal      - Column vector representing the complex or real signal.
%       sampleRate  - Sampling rate of the signal in Hz (positive scalar).
%
%   Outputs:
%       timeVector       - Time instants (seconds) corresponding to columns of W.
%       frequencyVector  - Non-negative frequency bins (Hz) corresponding to rows of W.
%       wvdMatrix        - Magnitude of the Wigner-Ville distribution.
%
%   Notes:
%       * The function rescales the lag FFT grid to keep only the
%         non-negative frequencies because the WVD of an analytic signal is
%         symmetric. This greatly reduces the memory footprint of the
%         distribution while preserving all energy content.
%       * The result is returned as the magnitude of the Wigner-Ville
%         distribution to facilitate downstream visualization and smoothing.
%
%   Example:
%       fs = 1e3;
%       t = (0:1023)'/fs;
%       x = exp(1j*2*pi*100*t);
%       [tVec, fVec, W] = computeWignerVilleDistribution(x, fs);
%
%   See also HILBERT, FFT.

    arguments
        signal (:, 1) double
        sampleRate (1, 1) double {mustBePositive}
    end

    % Ensure the signal is a column vector of doubles.
    x = double(signal(:));
    numSamples = length(x);

    if numSamples < 2
        error('computeWignerVilleDistribution:SignalTooShort', ...
              'The input signal must contain at least two samples.');
    end

    % Use the analytic signal to suppress negative-frequency aliasing.
    analyticSignal = hilbert(x);

    % Prepare the FFT grid in the lag dimension. Using a power-of-two grid
    % provides efficient FFT computation and balances time/frequency
    % resolution. Doubling the signal length ensures that all valid lags can
    % be represented without wrap-around artifacts.
    fftLength = 2^nextpow2(2 * numSamples);
    halfFftLength = fftLength / 2;

    lagFftBuffer = zeros(fftLength, 1);
    wvdSpectrum = zeros(fftLength, numSamples);

    % Compute the Wigner-Ville distribution by iterating over time samples
    % and evaluating the instantaneous autocorrelation sequence. For each
    % time instant we populate the lag-domain buffer symmetrically before
    % performing the FFT to obtain the frequency content.
    for n = 1:numSamples
        tauMax = min([n - 1, numSamples - n]);
        if tauMax == 0
            % Only the zero-lag term is valid at this time index.
            lagFftBuffer(:) = 0;
            lagFftBuffer(halfFftLength + 1) = abs(analyticSignal(n))^2;
        else
            tau = -tauMax:tauMax;
            lagValues = analyticSignal(n + tau) .* conj(analyticSignal(n - tau));

            lagFftBuffer(:) = 0;
            lagIndices = halfFftLength + 1 + tau;
            lagFftBuffer(lagIndices) = lagValues;
        end

        wvdSpectrum(:, n) = fft(lagFftBuffer);
    end

    % Align the spectrum so that zero frequency is centred and keep only the
    % real component (theoretically the WVD of an analytic signal is real).
    wvdSpectrum = real(fftshift(wvdSpectrum, 1));

    % Construct the associated time and frequency grids.
    fullFrequencyVector = linspace(-sampleRate / 2, sampleRate / 2, fftLength).';
    timeVector = (0:numSamples - 1).' / sampleRate;

    % Keep only the non-negative frequencies to reduce redundancy.
    positiveFrequencyMask = fullFrequencyVector >= 0;
    frequencyVector = fullFrequencyVector(positiveFrequencyMask) / 2;
    wvdMatrix = abs(wvdSpectrum(positiveFrequencyMask, :));
end


