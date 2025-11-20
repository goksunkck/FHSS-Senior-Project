% --- train_lstm_model.m ---
% 
% This script trains an LSTM network to detect FHSS hop sequences.
% It assumes that the data has already been prepared and partitioned by
% the 'prepare_dataset.m' script.
% 
% Workflow:
% 1. Load the pre-partitioned training and validation data.
% 2. Define the LSTM network architecture.
% 3. Specify training options.
% 4. Train the network.

clear; clc; close all; 

% --- 1. Load Prepared Dataset ---
fprintf('Loading prepared training and validation data...\n');
prepared_data_filename = 'data/synthetic/lstm_prepared_data.mat';
if ~exist(prepared_data_filename, 'file')
    error('Prepared dataset file not found. Please run prepare_dataset.m first.');
end
load(prepared_data_filename, 'XTrain', 'YTrain', 'XValidation', 'YValidation');

fprintf('Loaded %d training samples and %d validation samples.\n', ...
    length(XTrain), length(XValidation));

% --- 2. Define LSTM Network Architecture ---
fprintf('Defining LSTM network architecture...\n');

% Get feature dimension and number of classes from the loaded data.
if isempty(XTrain) || isempty(YTrain)
    error('Training data is empty. Cannot define network.');
end
featureDimension = size(XTrain{1}, 1); % Number of frequency bins in STFT
numClasses = numel(categories(YTrain{1})); % Number of unique hop indices

layers = [ ...
    sequenceInputLayer(featureDimension, 'Name', 'input')
    
    % A Bi-LSTM layer learns from both past and future context.
    % 'OutputMode','sequence' is crucial for sequence-to-sequence problems.
    bilstmLayer(128, 'OutputMode', 'sequence', 'Name', 'bilstm1')
    dropoutLayer(0.2, 'Name', 'dropout1')
    bilstmLayer(128, 'OutputMode', 'sequence', 'Name', 'bilstm2')
    dropoutLayer(0.2, 'Name', 'dropout2')

    % A fully connected layer processes each time step independently.
    fullyConnectedLayer(numClasses, 'Name', 'fc')
    
    % Softmax and classification layer for the final output.
    softmaxLayer('Name', 'softmax')
    classificationLayer('Name', 'classification')];

% You can analyze the network architecture using the Deep Network Designer app
% analyzeNetwork(layers);

% --- 3. Specify Training Options ---
fprintf('Specifying training options...\n');

% Set a mini-batch size. This is a key hyperparameter.
miniBatchSize = 32;

% Calculate how many iterations are in one epoch to validate once per epoch.
% Using length() is robust to row or column vectors of cells.
numIterationsPerEpoch = floor(length(XTrain) / miniBatchSize);
if numIterationsPerEpoch == 0
    numIterationsPerEpoch = 1; % Ensure it's at least 1
end

options = trainingOptions('adam', ...
    'ExecutionEnvironment', 'gpu', ...
    'MaxEpochs', 20, ...
    'MiniBatchSize', miniBatchSize, ...
    'ValidationData', {XValidation, YValidation}, ...
    'ValidationFrequency', numIterationsPerEpoch, ... % Validate once per epoch
    'GradientThreshold', 1, ...
    'InitialLearnRate', 0.001, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 5, ...
    'LearnRateDropFactor', 0.5, ...
    'SequenceLength', 'longest',... % Pad sequences to the length of the longest in the mini-batch
    'Shuffle', 'every-epoch', ...
    'Verbose', true, ...
    'Plots', 'training-progress');

% --- 4. Train the Network ---
fprintf('Script is ready for training.\n');
fprintf('This will open a new window showing the training progress.\n');
fprintf('Training may take a significant amount of time.\n');

[net, trainInfo] = trainNetwork(XTrain, YTrain, layers, options);

fprintf('\nTraining complete. Saving trained model...\n');
save('models/trained_lstm_model.mat', 'net', 'trainInfo');
fprintf('Trained model saved to models/trained_lstm_model.mat\n');

% After training, you can inspect the results:
% - 'net' contains the trained network.
% - 'trainInfo' contains information about the training process, like
%   final validation accuracy.
%
% Example: To classify a new sequence X_new:
% Y_pred = classify(net, X_new);
