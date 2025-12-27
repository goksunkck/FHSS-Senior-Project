% --- prepare_test_sequences.m ---
%
% Prepares sequences from the TEST dataset.
%
% Output: data/synthetic/prepared_test_sequences.mat

clear; clc; close all;

% --- Configuration ---
lookback_window = 35;

% --- Load TEST dataset ---
dataset_filename = 'data/synthetic/classification_dataset_stft_awgn_test.mat';
if ~exist(dataset_filename, 'file')
    error('Test dataset file not found. Run generate_test_dataset.m first.');
end

matObj = matfile(dataset_filename);
if isprop(matObj, 'numHopsPerSignal')
    numHopsPerSignal = matObj.numHopsPerSignal;
else
    numHopsPerSignal = 256;
end

num_total_examples = size(matObj, 'Y_data', 1);
num_signals = num_total_examples / numHopsPerSignal;

fprintf('Loaded %d total hops across %d signals.\n', num_total_examples, num_signals);

% --- Pre-allocate index + label arrays ---
num_sequences_per_signal = numHopsPerSignal - lookback_window;
total_sequences = num_signals * num_sequences_per_signal;
sequence_starts = zeros(total_sequences, 1, 'uint32');
Y_labels = zeros(total_sequences, 1, 'uint8');
sequence_counter = 0;

fprintf('Creating sliding windows (lookback = %d)...\n', lookback_window);
for i = 0:(num_signals - 1)
    signal_start_idx = i * numHopsPerSignal + 1;
    signal_end_idx = (i + 1) * numHopsPerSignal;

    current_Y_sequence = matObj.Y_data(signal_start_idx:signal_end_idx, 1);

    for j = 1:num_sequences_per_signal
        sequence_counter = sequence_counter + 1;
        % sequence_start is the index of the first hop in the lookback window
        sequence_starts(sequence_counter) = signal_start_idx + j - 1;
        % label is the hop immediately AFTER the window
        Y_labels(sequence_counter) = current_Y_sequence(j + lookback_window);
    end
end

fprintf('Generated %d sequences.\n', sequence_counter);

% --- Metadata ---
k = 3;
numChannels = 2^k;
class_values = uint8(0:numChannels-1)';
class_names = string(class_values);

% No train/val split for test set, we use all valid sequences.
idxValidation = true(total_sequences, 1);
sequence_starts_validation = sequence_starts;
YValidation = Y_labels;

% --- Save prepared dataset ---
output_filename = 'data/synthetic/prepared_test_sequences.mat';
fprintf('Saving prepared TEST sequences to %s...\n', output_filename);
save(output_filename, 'sequence_starts', 'sequence_starts_validation', ...
    'Y_labels', 'YValidation', 'idxValidation', ...
    'class_values', 'class_names', 'lookback_window', 'numHopsPerSignal', ...
    'dataset_filename', '-v7.3');

fprintf('Preparation complete.\n');
