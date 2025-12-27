% --- prepare_prediction_sequences.m ---
%
% Converts the per-hop STFT dataset into an index-based forecasting dataset.
% Stores only:
%   - sequence_starts: start indices of lookback windows
%   - labels: next-hop indices
%   - class mapping for categorical casting
% The STFT tensors remain in the source MAT-file and are loaded on-demand.
%
% Output: data/synthetic/prepared_prediction_sequences.mat

clear; clc; close all;
% TODO: Make sure to modify prepare_prediction_sequences.m to split the data into training: first 8 hopsets, validation: last 2 hopsets
% --- Configuration ---
lookback_window = 35;      % Number of previous hops to use for prediction
validation_split_ratio = 0.2;

% --- Load base dataset ---
dataset_filename = 'data/synthetic/classification_dataset_stft_awgn.mat';
if ~exist(dataset_filename, 'file')
    error('Dataset file not found. Run generate_prediction_dataset.m first.');
end

matObj = matfile(dataset_filename);
if isprop(matObj, 'numHopsPerSignal')
    numHopsPerSignal = matObj.numHopsPerSignal;
else
    numHopsPerSignal = 256; % fallback if absent
end

num_total_examples = size(matObj, 'Y_data', 1);
num_signals = num_total_examples / numHopsPerSignal;
if mod(num_signals, 1) ~= 0
    error('Total examples (%d) not divisible by numHopsPerSignal (%d). Check dataset.', ...
        num_total_examples, numHopsPerSignal);
end
if lookback_window >= numHopsPerSignal
    error('lookback_window (%d) must be smaller than numHopsPerSignal (%d).', ...
        lookback_window, numHopsPerSignal);
end

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
        sequence_starts(sequence_counter) = signal_start_idx + j - 1;
        Y_labels(sequence_counter) = current_Y_sequence(j + lookback_window);
    end

    if mod(i+1, 20) == 0
        fprintf('  Processed %d / %d signals...\n', i+1, num_signals);
    end
end

fprintf('Generated %d sequences.\n', sequence_counter);

% --- Prepare labels and split ---
% --- Prepare labels and split ---
% Define classes based on system parameters
k = 3;
numChannels = 2^k;
class_values = uint8(0:numChannels-1)';
class_names = string(class_values);

% --- Time Split Strategy (Tracking) ---
% We use the first 80% of EVERY signal for training (Learning the Pattern)
% We use the last 20% of EVERY signal for validation (Tracking the Pattern)
split_idx = floor(num_sequences_per_signal * 0.8);

fprintf('Splitting data by TIME (Tracking Mode)...\n');
fprintf('  Training: First %d sequences of each signal.\n', split_idx);
fprintf('  Validation: Remaining %d sequences of each signal.\n', num_sequences_per_signal - split_idx);

idxTrain_list = false(total_sequences, 1);
idxValidation_list = false(total_sequences, 1);

current_seq_idx = 0;
for i = 0:(num_signals - 1)
    for j = 1:num_sequences_per_signal
        current_seq_idx = current_seq_idx + 1;

        % If j <= split_idx -> Train
        if j <= split_idx
            idxTrain_list(current_seq_idx) = true;
        else
            idxValidation_list(current_seq_idx) = true;
        end
    end
end

idxTrain = idxTrain_list;
idxValidation = idxValidation_list;

sequence_starts_train = sequence_starts(idxTrain);
YTrain = Y_labels(idxTrain);
sequence_starts_validation = sequence_starts(idxValidation);
YValidation = Y_labels(idxValidation);

fprintf('Time Split Complete.\n');
fprintf('  Train sequences: %d (%.1f%%)\n', numel(YTrain), 100*numel(YTrain)/total_sequences);
fprintf('  Validation sequences: %d (%.1f%%)\n', numel(YValidation), 100*numel(YValidation)/total_sequences);

% --- Save prepared dataset ---
output_filename = 'data/synthetic/prepared_prediction_sequences.mat';
fprintf('Saving prepared prediction sequences to %s...\n', output_filename);
save(output_filename, 'sequence_starts', 'sequence_starts_train', 'sequence_starts_validation', ...
    'Y_labels', 'YTrain', 'YValidation', 'idxTrain', 'idxValidation', ...
    'class_values', 'class_names', 'lookback_window', 'numHopsPerSignal', ...
    'dataset_filename', '-v7.3');

fprintf('Preparation complete.\n');
