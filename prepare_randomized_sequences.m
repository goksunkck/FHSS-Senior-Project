% --- prepare_randomized_sequences.m ---
%
% Prepares sequences from the RANDOMIZED DEGREE 10 dataset.
% Implements a Signal-wise Split (70% Train, 15% Val, 15% Test).
%
% Output: data/synthetic/prepared_randomized_sequences.mat

clear; clc; close all;

% --- Configuration ---
lookback_window = 35;
dataset_filename = 'data/synthetic/classification_dataset_stft_random_deg10.mat';

if ~exist(dataset_filename, 'file')
    error('Dataset file not found. Run generate_randomized_dataset.m first.');
end

matObj = matfile(dataset_filename);
numHopsPerSignal = matObj.numHopsPerSignal;
num_total_examples = size(matObj, 'Y_data', 1);
num_signals = num_total_examples / numHopsPerSignal;

fprintf('Loaded %d total hops across %d signals.\n', num_total_examples, num_signals);

% --- Pre-allocate ---
num_sequences_per_signal = numHopsPerSignal - lookback_window;
total_sequences = num_signals * num_sequences_per_signal;
sequence_starts = zeros(total_sequences, 1, 'uint32');
Y_labels = zeros(total_sequences, 1, 'uint8');

% --- Track Signal IDs for splitting ---
% We store which signal each sequence belongs to
signal_ids = zeros(total_sequences, 1, 'uint32');

sequence_counter = 0;

fprintf('Creating sliding windows...\n');
for i = 0:(num_signals - 1)
    signal_start_idx = i * numHopsPerSignal + 1;
    signal_end_idx = (i + 1) * numHopsPerSignal;

    % Load signal's Label Sequence
    current_Y_sequence = matObj.Y_data(signal_start_idx:signal_end_idx, 1);

    for j = 1:num_sequences_per_signal
        sequence_counter = sequence_counter + 1;
        sequence_starts(sequence_counter) = signal_start_idx + j - 1;
        Y_labels(sequence_counter) = current_Y_sequence(j + lookback_window);
        signal_ids(sequence_counter) = i; % Signal ID (0-indexed)
    end

    if mod(i, 100) == 0
        fprintf('  Processed %d signals...\n', i);
    end
end

% --- Splitting Logic (Signal-wise) ---
% We split signals, NOT sequences. This ensures no "sequence leakage"
% from a training signal into validation.
all_signal_indices = unique(signal_ids);
num_unique_signals = length(all_signal_indices);

% Random shuffle of signal IDs
rng(42); % Fixed seed for reproducibility of split
shuffled_signals = all_signal_indices(randperm(num_unique_signals));

% 70% Train, 15% Val, 15% Test
num_train = floor(0.7 * num_unique_signals);
num_val = floor(0.15 * num_unique_signals);
num_test = num_unique_signals - num_train - num_val;

train_signals = shuffled_signals(1:num_train);
val_signals = shuffled_signals(num_train+1 : num_train+num_val);
test_signals = shuffled_signals(num_train+num_val+1 : end);

fprintf('Split: Train=%d, Val=%d, Test=%d signals.\n', length(train_signals), length(val_signals), length(test_signals));

% Create Boolean Masks
% Ideally, we map signal_id -> set logic.
% Since signal_ids are integers 0..N, we can use a lookup table.
split_map = zeros(num_unique_signals, 1, 'uint8'); % 0:Train, 1:Val, 2:Test
split_map(val_signals + 1) = 1; % indices are 1-based in MATLAB
split_map(test_signals + 1) = 2;

% Map to sequences
sequence_splits = split_map(signal_ids + 1);

idxTrain = sequence_splits == 0;
idxValidation = sequence_splits == 1;
idxTest = sequence_splits == 2;

% --- Labels ---
k = 3;
numChannels = 2^k;
class_values = uint8(0:numChannels-1)';
class_names = string(class_values);

% --- Select Subsets ---
% We construct the final arrays to save
sequence_starts_train = sequence_starts(idxTrain);
YTrain = Y_labels(idxTrain);

sequence_starts_validation = sequence_starts(idxValidation);
YValidation = Y_labels(idxValidation);

sequence_starts_test = sequence_starts(idxTest);
YTest = Y_labels(idxTest);

% --- Global Stats ---
% Copy global stats from source dataset to prepared file
global_min = -100;
global_max = 100;
try
    global_min = matObj.global_min;
    global_max = matObj.global_max;
catch
    warning('Global stats not found in source dataset.');
end

% --- Save ---
output_filename = 'data/synthetic/prepared_randomized_sequences.mat';
fprintf('Saving prepared sequences to %s...\n', output_filename);
save(output_filename, 'sequence_starts', 'sequence_starts_train', 'sequence_starts_validation', 'sequence_starts_test', ...
    'Y_labels', 'YTrain', 'YValidation', 'YTest', ...
    'idxTrain', 'idxValidation', 'idxTest', ...
    'class_values', 'class_names', 'lookback_window', 'numHopsPerSignal', ...
    'dataset_filename', 'global_min', 'global_max', '-v7.3');

fprintf('Preparation complete.\n');
