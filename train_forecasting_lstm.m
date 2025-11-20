% --- train_forecasting_lstm.m ---
% 
% This script trains an LSTM network to FORECAST the next hop in a
% frequency-hopping sequence, based on a window of previous hops.
% 
% Workflow:
% 1. Load the pre-partitioned forecasting data.
% 2. Define the LSTM network architecture for forecasting.
% 3. Specify training options.
% 4. Train the network.
% 5. Save the trained model.

clear; clc; close all;

% --- 1. Load Prepared Forecasting Dataset ---
fprintf('Loading prepared forecasting data...\n');
prepared_data_filename = 'data/synthetic/lstm_forecasting_data.mat';
if ~exist(prepared_data_filename, 'file')
    error('Prepared dataset file not found. Please run prepare_forecasting_dataset.m first.');
end
load(prepared_data_filename, 'XTrain', 'YTrain', 'XValidation', 'YValidation');

fprintf('Loaded %d training samples and %d validation samples.\n', ...
    length(XTrain), length(XValidation));

% --- 2. Define LSTM Network Architecture for Forecasting ---
fprintf('Defining LSTM network architecture for forecasting...\n');

if isempty(YTrain)
    error('Training data is empty. Cannot define network.');
end
numClasses = numel(categories(YTrain)); % Number of unique hop indices

% NOTE: The input dimension is 1 because we are feeding in one hop index
% at a time from the lookback window.
layers = [ ...
    sequenceInputLayer(1, 'Name', 'input')

    % Use standard LSTM layers for forecasting, not Bi-LSTM.
    % 'OutputMode','sequence' passes the full sequence to the next layer.
    lstmLayer(128, 'OutputMode', 'sequence', 'Name', 'lstm1')
    dropoutLayer(0.2, 'Name', 'dropout1')
    
    % The second LSTM layer outputs only the *last* time step, which is
    % what we need for a sequence-to-one prediction.
    lstmLayer(128, 'OutputMode', 'last', 'Name', 'lstm2')
    dropoutLayer(0.2, 'Name', 'dropout2')

    % A fully connected layer processes the output of the last time step.
    fullyConnectedLayer(numClasses, 'Name', 'fc')
    
    % Softmax and classification layer for the final output.
    softmaxLayer('Name', 'softmax')
    classificationLayer('Name', 'classification')];

% You can analyze the network architecture using the Deep Network Designer app
% analyzeNetwork(layers);

% --- 3. Specify Training Options ---
fprintf('Specifying training options...\n');

options = trainingOptions('adam', ...
    'ExecutionEnvironment', 'gpu', ... % Explicitly use 'gpu'
    'MaxEpochs', 30, ... % May need more epochs for this kind of task
    'MiniBatchSize', 128, ...
    'ValidationData', {XValidation, YValidation}, ...
    'ValidationFrequency', 500, ...
    'GradientThreshold', 1, ...
    'InitialLearnRate', 0.002, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 10, ...
    'LearnRateDropFactor', 0.5, ...
    'Shuffle', 'every-epoch', ...
    'Verbose', true, ...
    'Plots', 'training-progress');

% --- 4. Train the Network ---
fprintf('Script is ready for training the forecasting model.\n');
fprintf('This will open a new window showing the training progress.\n');

[net, trainInfo] = trainNetwork(XTrain, YTrain, layers, options);

% --- 5. Save the Trained Model ---
fprintf('\nTraining complete. Saving trained forecasting model...\n');
model_output_filename = 'models/forecasting_lstm_model.mat';
save(model_output_filename, 'net', 'trainInfo');
fprintf('Trained forecasting model saved to %s\n', model_output_filename);

% Example: To predict the next hop for a new sequence X_new (where X_new
% is a column vector of previous hops):
% next_hop_pred = classify(net, {X_new});
