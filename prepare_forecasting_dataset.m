% --- prepare_forecasting_dataset.m ---
% 
% This script transforms the sequence of hop indices into a dataset
% suitable for time-series forecasting. It uses a sliding window approach
% to create input-output pairs.
% 
% Workflow:
% 1. Load the raw dataset containing the hop sequences (Y_data).
% 2. Define a 'lookback_window' size.
% 3. For each sequence, create pairs of [sequence_window, next_hop].
% 4. Partition the new data into training and validation sets.
% 5. Save the prepared data to a new .mat file.
% 
% Output file: 'data/synthetic/lstm_forecasting_data.mat'

clear; clc; close all;

% --- 1. Configuration ---
lookback_window = 15; % Number of previous hops to use for prediction
validation_split_ratio = 0.2; % 20% of data for validation

% --- 2. Load Raw Hop Sequences ---
fprintf('Loading raw dataset for hop sequences...\n');
dataset_filename = 'data/synthetic/lstm_dataset_stft_awgn.mat';
if ~exist(dataset_filename, 'file')
    error('Raw dataset file not found. Please run create_lstm_dataset_stft.m first.');
end
% We only need Y_data for this task.
load(dataset_filename, 'Y_data');

fprintf('Loaded %d raw hop sequences.\n', length(Y_data));

% --- 3. Create Sliding Window Dataset ---
fprintf('Creating sliding window dataset with a lookback of %d...\n', lookback_window);

X_sequences = {}; % Cell array for input sequences
Y_labels = [];    % Vector for corresponding labels

for i = 1:length(Y_data)
    sequence = Y_data{i};
    % Slide a window across each sequence
    for j = 1:(length(sequence) - lookback_window)
        
        % The input is a window of N hops
        input_window = sequence(j : j + lookback_window - 1);
        
        % The label is the very next hop
        output_label = sequence(j + lookback_window);
        
        % Store the pair
        X_sequences{end+1} = input_window'; % Store as column vector
        Y_labels(end+1) = output_label;
        
    end
end

fprintf('Generated %d training examples.\n', length(X_sequences));

% --- 4. Prepare Labels and Partition Data ---

% Convert labels to a categorical array. The +1 is because our data is
% 0-based but MATLAB categoricals are 1-based.
numClasses = max(Y_labels) + 1;
classNames = string(0:(numClasses-1));
Y_labels_cat = categorical(Y_labels, 0:(numClasses-1), classNames);

fprintf('Partitioning data into training and validation sets (%.0f/%.0f split)...\n', ...
    (1-validation_split_ratio)*100, validation_split_ratio*100);

% Set the random seed for reproducibility
rng(42);

% Create a cross-validation partition object.
cv = cvpartition(length(X_sequences), 'HoldOut', validation_split_ratio);
idxTrain = training(cv);
idxValidation = test(cv);

% Select the data for each set.
XTrain = X_sequences(idxTrain);
YTrain = Y_labels_cat(idxTrain);
XValidation = X_sequences(idxValidation);
YValidation = Y_labels_cat(idxValidation);

% Transpose YTrain to be a column vector of categoricals, as required by trainNetwork
YTrain = YTrain';
YValidation = YValidation';

fprintf('Data partitioned:\n');
fprintf('  - Training samples:   %d\n', length(XTrain));
fprintf('  - Validation samples: %d\n', length(XValidation));


% --- 5. Save Prepared Data ---
output_filename = 'data/synthetic/lstm_forecasting_data.mat';
fprintf('Saving prepared forecasting data to %s...\n', output_filename);

save(output_filename, 'XTrain', 'YTrain', 'XValidation', 'YValidation', '-v7.3');

fprintf('Data preparation for forecasting complete.\n');
