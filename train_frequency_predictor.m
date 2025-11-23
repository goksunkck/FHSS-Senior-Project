% --- train_frequency_predictor.m ---
%
% Trains a CNN-LSTM to predict the next hop frequency index using
% in-memory loading of STFT windows from the base dataset.
%
% Requires:
%   - data/synthetic/prepared_prediction_sequences.mat (index + label sets)
%   - data/synthetic/classification_dataset_stft_awgn.mat (STFT tensors)
%
% Output:
%   - models/trained_frequency_predictor.mat

clear; clc; close all;

% --- 1. Load prepared data ---
% This script assumes 'prepare_prediction_sequences.m' has been run with the
% correct, robust class definitions and has produced a clean data file.
prep_file = 'data/synthetic/prepared_prediction_sequences.mat';
if ~exist(prep_file, 'file')
    error('Prepared dataset not found. Please delete it and re-run prepare_prediction_sequences.m first.');
end
load(prep_file, 'sequence_starts_train', 'sequence_starts_validation', ...
    'YTrain', 'YValidation', 'lookback_window', 'dataset_filename', ...
    'class_values', 'class_names');

% --- 2. Load base STFT dataset into memory ---
source_dataset = dataset_filename;
if ~exist(source_dataset, 'file')
    error('Source dataset not found at %s. Please re-run generate_prediction_dataset.m.', source_dataset);
end
fprintf('Loading entire STFT dataset (%s) into memory...\n', source_dataset);
data_contents = load(source_dataset, 'X_data');
X_data_in_memory = data_contents.X_data;
clear data_contents;
fprintf('Dataset loaded into memory.\n');

fprintf('Loaded %d training and %d validation sequence indices.\n', ...
    numel(sequence_starts_train), numel(sequence_starts_validation));

% --- 3. Process Labels ---
% Define the master list of categories from the clean file
all_categories = string(class_names);

% Convert numeric labels from the file into categorical arrays
YTrain = categorical(YTrain, class_values, all_categories);
YValidation = categorical(YValidation, class_values, all_categories);

% Final integrity check. If this fails, the .mat file is corrupt.
if any(isundefined(YTrain)) || any(isundefined(YValidation))
    error('Undefined labels detected. Please delete prepared_prediction_sequences.mat and re-run the script to generate it.');
end

% --- 4. Determine input size ---
if isempty(sequence_starts_train)
    error('The training dataset is empty. Cannot proceed.');
end
firstCell = X_data_in_memory(sequence_starts_train(1), 1);
firstImage = firstCell{1};
firstSize = size(firstImage);
if numel(firstSize) < 3, firstSize(3) = 1; end
inputSize = [firstSize(1), firstSize(2), firstSize(3)];
numClasses = numel(all_categories);

fprintf('Input STFT size: [%d %d %d]\n', inputSize(1), inputSize(2), inputSize(3));
fprintf('Num classes: %d\n', numClasses);

% --- 5. Define CNN-LSTM architecture ---
layers = [
    sequenceInputLayer(inputSize, 'Name', 'input')
    convolution2dLayer(3, 16, 'Padding', 'same', 'Name', 'conv1')
    batchNormalizationLayer('Name', 'bn1')
    reluLayer('Name', 'relu1')
    maxPooling2dLayer([2 1], 'Stride', [2 1], 'Name', 'pool1')
    convolution2dLayer(3, 32, 'Padding', 'same', 'Name', 'conv2')
    batchNormalizationLayer('Name', 'bn2')
    reluLayer('Name', 'relu2')
    maxPooling2dLayer([2 1], 'Stride', [2 1], 'Name', 'pool2')
    flattenLayer('Name', 'flatten')
    lstmLayer(100, 'OutputMode', 'last', 'Name', 'lstm')
    fullyConnectedLayer(numClasses, 'Name', 'fc')
    softmaxLayer('Name', 'softmax')
    classificationLayer('Name', 'classification')
];

% --- 6. Datastores (in-memory sequence loading) ---
trainDs = createSequenceDatastore(X_data_in_memory, sequence_starts_train, YTrain, lookback_window, inputSize);
valDs = createSequenceDatastore(X_data_in_memory, sequence_starts_validation, YValidation, lookback_window, inputSize);

% --- 7. Training options ---
options = trainingOptions('adam', ...
    'ExecutionEnvironment', 'cpu', ...
    'InitialLearnRate', 1e-3, ...
    'MaxEpochs', 5, ...
    'MiniBatchSize', 4, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', valDs, ...
    'ValidationFrequency', 100, ...
    'DispatchInBackground', false, ...
    'Verbose', true, ...
    'Plots', 'training-progress');

% --- 8. Train ---
fprintf('Starting training...\n');
[trainedNet, trainInfo] = trainNetwork(trainDs, layers, options);
fprintf('Training complete.\n');

% --- 9. Save model ---
output_filename = 'models/trained_frequency_predictor.mat';
fprintf('Saving trained model to %s...\n', output_filename);
save(output_filename, 'trainedNet', 'trainInfo', 'inputSize', 'numClasses', ...
    'class_values', 'class_names', 'lookback_window');

fprintf('All done.\n');

% --- Local helper functions ---
function ds = createSequenceDatastore(X_data, sequenceStarts, labels, lookbackWindow, inputSize)
    % This helper creates a datastore that combines sequence start indices
    % and labels, and then transforms it to load the actual image data.
    ds_in = combine(arrayDatastore(sequenceStarts), arrayDatastore(labels));
    ds = transform(ds_in, @(data) localLoadSequence(data, X_data, lookbackWindow, inputSize));
end

function [X, Y] = localLoadSequence(data, X_data, lookbackWindow, inputSize)
    % This function processes one sample at a time. `trainNetwork` will
    % automatically batch the output of this function.
    % 'data' is a 1x2 cell {startIndex, label}.
    startIdx = data{1};
    Y = data{2}; % Label is already the correct categorical type

    idxRange = startIdx : (startIdx + lookbackWindow - 1);
    
    seqMatrix = zeros([inputSize, lookbackWindow], 'single');
    for t = 1:lookbackWindow
        imgCell = X_data(idxRange(t), 1);
        img = imgCell{1};
        
        % Brute-force reshape: Unconditionally reshape the image to be 3D.
        % This is safe because we know the source data is 2D.
        img = reshape(img, [size(img, 1), size(img, 2), 1]);
        
        % Validate the final size
        if ~isequal(size(img), inputSize)
             error('Input size mismatch at index %d after reshape. Expected [%s], but got size [%s].', ...
                idxRange(t), mat2str(inputSize), mat2str(size(img)));
        end
        seqMatrix(:, :, :, t) = img;
    end
    
    X = seqMatrix;
end

