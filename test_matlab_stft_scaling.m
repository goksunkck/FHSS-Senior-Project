% test_matlab_stft_scaling.m
% Quick test to see MATLAB's spectrogram scaling

fs = 10e6;
duration = 0.001; % 1 ms
t = (0:1/fs:duration-1/fs)';

% Create unit amplitude complex tone at 1 MHz
freq = 1e6;
signal = exp(1j * 2 * pi * freq * t);

fprintf('Signal length: %d\n', length(signal));
fprintf('Signal power: %.6f\n', mean(abs(signal).^2));
fprintf('Signal power (dB): %.2f\n', 10*log10(mean(abs(signal).^2)));

% STFT parameters
window_len = 128;
overlap_len = 96;
nfft = 1024;

% Compute spectrogram
[S, F, T] = spectrogram(signal, hamming(window_len), overlap_len, nfft, fs, 'centered');

% Log magnitude
mag = abs(S);
log_mag = 20*log10(mag + eps);

fprintf('\nMATLAB Spectrogram Output:\n');
fprintf('  Shape: [%d x %d]\n', size(S,1), size(S,2));
fprintf('  Max magnitude (linear): %.6f\n', max(mag(:)));
fprintf('  Max magnitude (dB): %.2f\n', max(log_mag(:)));
fprintf('  Mean magnitude (dB): %.2f\n', mean(log_mag(:)));

fprintf('\n========================================\n');
fprintf('EXPECTED vs ACTUAL:\n');
fprintf('========================================\n');
fprintf('Input signal power: 0 dB (unit amplitude)\n');
fprintf('Spectrogram peak magnitude: %.2f dB\n', max(log_mag(:)));
fprintf('Difference: %.2f dB\n', max(log_mag(:)) - 0);
