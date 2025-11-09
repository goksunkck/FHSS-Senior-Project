function [timeVector, frequencyVector, spwvdMatrix] = computeSmoothedPseudoWignerVille(signal, sampleRate, varargin)
%COMPUTESMOOTHEDPSEUDOWIGNERVILLE Smoothed pseudo Wigner-Ville distribution.
%   [T, F, S] = COMPUTESMOOTHEDPSEUDOWIGNERVILLE(X, FS) returns the time
%   vector T, frequency vector F, and smoothed pseudo Wigner-Ville matrix S
%   for the input signal X sampled at FS Hz. The function can accept a
%   pre-computed Wigner-Ville distribution to avoid redundant computation
%   when both the unsmoothed and smoothed representations are required.
%
%   Additional name-value pairs:
%       'TimeWindow'       - Vector defining the time smoothing window.
%                            Defaults to a length-31 Hann window.
%       'FrequencyWindow'  - Vector defining the frequency smoothing window.
%                            Defaults to a length-31 Hann window.
%       'WVDMatrix'        - Precomputed Wigner-Ville distribution. If
%                            provided, 'TimeVector' and 'FrequencyVector'
%                            must also be specified.
%       'TimeVector'       - Time instants corresponding to columns of the
%                            supplied Wigner-Ville distribution.
%       'FrequencyVector'  - Frequency bins corresponding to rows of the
%                            supplied Wigner-Ville distribution.
%
%   Example:
%       [tVec, fVec, W] = computeWignerVilleDistribution(x, fs);
%       [tVec, fVec, S] = computeSmoothedPseudoWignerVille([], fs, ...
%                           'WVDMatrix', W, 'TimeVector', tVec, ...
%                           'FrequencyVector', fVec);
%
%   See also COMPUTEWIGNERVILLEDISTRIBUTION, CONV2, HANN.

    parser = inputParser;
    parser.addParameter('TimeWindow', hann(31), @(w) isnumeric(w) && isvector(w));
    parser.addParameter('FrequencyWindow', hann(31), @(w) isnumeric(w) && isvector(w));
    parser.addParameter('WVDMatrix', []);
    parser.addParameter('TimeVector', []);
    parser.addParameter('FrequencyVector', []);
    parser.parse(varargin{:});

    timeWindow = parser.Results.TimeWindow(:);
    freqWindow = parser.Results.FrequencyWindow(:);
    wvdMatrix = parser.Results.WVDMatrix;
    timeVector = parser.Results.TimeVector;
    frequencyVector = parser.Results.FrequencyVector;

    if isempty(wvdMatrix)
        [timeVector, frequencyVector, wvdMatrix] = computeWignerVilleDistribution(signal, sampleRate);
    else
        if isempty(timeVector) || isempty(frequencyVector)
            error('computeSmoothedPseudoWignerVille:MissingAxes', ...
                  ['TimeVector and FrequencyVector must be provided when ', ...
                   'supplying a precomputed Wigner-Ville distribution.']);
        end
    end

    % Normalise the smoothing windows so that the overall energy is
    % preserved. Guard against zero-valued windows.
    timeWindow = normaliseWindow(timeWindow, 'TimeWindow');
    freqWindow = normaliseWindow(freqWindow, 'FrequencyWindow');

    % Construct the two-dimensional smoothing kernel via the outer product
    % of the individual time and frequency windows.
    smoothingKernel = freqWindow * timeWindow.';

    % Apply the smoothing kernel to the Wigner-Ville distribution. The
    % result is guaranteed to remain real-valued because both the kernel and
    % the original distribution are real.
    spwvdMatrix = conv2(wvdMatrix, smoothingKernel, 'same');
    spwvdMatrix = max(spwvdMatrix, 0);
end

function window = normaliseWindow(window, name)
%NORMALISEWINDOW Normalise a window function and validate its contents.
    if isempty(window)
        error('computeSmoothedPseudoWignerVille:EmptyWindow', ...
              '%s must be a non-empty vector.', name);
    end

    window = window(:);
    totalEnergy = sum(window);
    if totalEnergy == 0
        error('computeSmoothedPseudoWignerVille:ZeroWindow', ...
              '%s cannot be identically zero.', name);
    end

    window = window / totalEnergy;
end