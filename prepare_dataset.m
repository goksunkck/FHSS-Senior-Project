% --- prepare_dataset.m ---
% 
% This script loads the generated FHSS dataset, converts labels to a
% categorical format, partitions it into training and validation sets,
% and saves the prepared data to a new .mat file.
% 
% Output file: 'data/synthetic/lstm_prepared_data.mat' which contains:
%   - XTrain:      Cell array of STFT matrices for training.
%   - YTrain:      Cell array of categorical label sequences for training.
%   - XValidation: Cell array of STFT matrices for validation.
%   - YValidation: Cell array of categorical label sequences for validation.

clear; clc; close all; 

% --- 1. Load Raw Dataset ---
fprintf('Loading raw dataset...\n');
dataset_filename = 'data/synthetic/lstm_dataset_stft_awgn.mat';
if ~exist(dataset_filename, 'file')
    error('Dataset file not found. Please run create_lstm_dataset_stft.m first.');
end
load(dataset_filename, 'X_data', 'Y_data');

fprintf('Loaded %d raw signals.\n', length(X_data));

% --- 2. Prepare Labels for Classification ---
fprintf('Aligning labels to spectrogram time steps for sequence-to-sequence learning...\n');

% These parameters from create_lstm_dataset_stft.m are needed for alignment.
stft_params.windowLength = 256;
stft_params.overlapLength = 192;
simParams.numModulatedSamples = 51200; % From 512 bits, 100 samples/symbol in the generator
simParams.numHops = 256;

% Calculate STFT hop size (time resolution) and samples per signal hop
stft_hop_size = stft_params.windowLength - stft_params.overlapLength; % Should be 64
samples_per_hop = floor(simParams.numModulatedSamples / simParams.numHops); % Should be 200

if isempty(Y_data)
    error('Y_data is empty. No data to process.');
end

Y_data_aligned = cell(size(Y_data));
for i = 1:length(X_data)
    spectrogram_seq_len = size(X_data{i}, 2); % Number of columns (time steps)
    original_hop_sequence = Y_data{i}; % The original 256 hops
    
    aligned_labels = zeros(1, spectrogram_seq_len);
    for t = 1:spectrogram_seq_len
        % Calculate the sample index at the center of this STFT frame
        center_sample_index = (t - 1) * stft_hop_size + floor(stft_params.windowLength / 2);
        
        % Determine which signal hop this sample belongs to
        hop_index_for_this_frame = floor((center_sample_index - 1) / samples_per_hop) + 1;
        
        % Clamp the index to be within the bounds of the original hop sequence
        hop_index_for_this_frame = max(1, min(hop_index_for_this_frame, simParams.numHops));
        
        % Assign the true hop value from the original sequence
        aligned_labels(t) = original_hop_sequence(hop_index_for_this_frame);
    end
    Y_data_aligned{i} = aligned_labels;
end

% Now, convert the aligned numerical sequences to categorical sequences
% The +1 is because MATLAB categories are 1-based, but our data is 0-based.
numClasses = max(cellfun(@max, Y_data_aligned)) + 1;
classNames = string(0:(numClasses-1));

Y_data_cat = cell(size(Y_data_aligned));
for i = 1:length(Y_data_aligned)
    % The third argument 'classNames' ensures all categoricals have the same set of classes
    Y_data_cat{i} = categorical(Y_data_aligned{i}, 0:(numClasses-1), classNames);
end

% --- 3. Partition Data ---
fprintf('Partitioning data into training and validation sets (80/20 split)...\n');

% Set the random seed for reproducibility of the split
rng(42); 

% Create a cross-validation partition object.
cv = cvpartition(length(X_data), 'HoldOut', 0.2);
idxTrain = training(cv);
idxValidation = test(cv);

% Select the data for each set based on the indices.
XTrain = X_data(idxTrain);
YTrain = Y_data_cat(idxTrain);
XValidation = X_data(idxValidation);
YValidation = Y_data_cat(idxValidation);

fprintf('Data partitioned:\n');
fprintf('  - Training samples:   %d\n', length(XTrain));
fprintf('  - Validation samples: %d\n', length(XValidation));

% --- 4. Save Prepared Data ---
output_filename = 'data/synthetic/lstm_prepared_data.mat';
fprintf('Saving prepared data to %s...\n', output_filename);

save(output_filename, 'XTrain', 'YTrain', 'XValidation', 'YValidation', '-v7.3');

fprintf('Data preparation complete.\n');
