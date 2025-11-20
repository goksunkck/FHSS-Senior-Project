% --- test_lstm_model_by_snr.m ---
%
% This script evaluates the trained LSTM model's performance at various SNR levels.
%
% Workflow:
% 1. Load the trained LSTM model.
% 2. Load the complete, raw dataset containing signals for all SNRs.
% 3. Reconstruct which signals belong to which SNR based on the generation script's logic.
% 4. For each SNR level:
%    a. Isolate the test data (X and Y) for that SNR.
%    b. Pre-process the labels (Y) to be categorical.
%    c. Use the trained network to classify the signals.
%    d. Calculate the classification accuracy.
% 5. Plot the final accuracy vs. SNR results.

clear; clc; close all;
addpath(genpath('src'));

% --- 1. Load Trained Model ---
fprintf('Loading trained LSTM model...\n');
model_filename = 'models/trained_lstm_model.mat';
if ~exist(model_filename, 'file')
    error('Trained model file not found. Please run train_lstm_model.m first.');
end
load(model_filename, 'net');

% --- 2. Load Raw Dataset ---
fprintf('Loading raw dataset...\n');
dataset_filename = 'data/synthetic/lstm_dataset_stft_awgn.mat';
if ~exist(dataset_filename, 'file')
    error('Raw dataset file not found. Please run create_lstm_dataset_stft.m first.');
end
load(dataset_filename, 'X_data', 'Y_data');
fprintf('Loaded %d total signals.\n', length(X_data));

% --- 3. Reconstruct SNR Information ---
% These parameters MUST match the ones in create_lstm_dataset_stft.m
snr_levels_db = -10:2:10;
num_hopsets = 10;
num_signals_per_config = 9;
signals_per_snr = num_hopsets * num_signals_per_config;

if length(X_data) ~= length(snr_levels_db) * signals_per_snr
    error('The number of signals in the dataset does not match the expected number based on the generation script parameters. Please check the parameters in both scripts.');
end

% Define hopset to map indices to frequencies (must match generation script)
numChannels = 8; % 2^k where k=3
freqSeparation = 10e6 * 0.01; % fs * 0.01
spacing = 2 * freqSeparation;
baseFreq = 2e6;
hopset = (0:numChannels-1) * spacing + baseFreq;


% --- 4. Evaluate Performance per SNR ---
fprintf('Starting evaluation for each SNR...\n');
accuracies = zeros(size(snr_levels_db));

% Label alignment parameters (must match prepare_dataset.m)
stft_params.windowLength = 256;
stft_params.overlapLength = 192;
simParams.numModulatedSamples = 51200;
simParams.numHops = 256;
stft_hop_size = stft_params.windowLength - stft_params.overlapLength;
samples_per_hop = floor(simParams.numModulatedSamples / simParams.numHops);

for i = 1:length(snr_levels_db)
    snr_db = snr_levels_db(i);
    fprintf('--- Evaluating SNR = %d dB ---\n', snr_db);

    % a. Isolate data for the current SNR
    start_idx = (i-1) * signals_per_snr + 1;
    end_idx = i * signals_per_snr;
    X_test_snr = X_data(start_idx:end_idx);
    Y_test_snr_raw = Y_data(start_idx:end_idx);

    % b. Pre-process labels for this SNR's data
    Y_test_aligned = cell(size(Y_test_snr_raw));
    for j = 1:length(X_test_snr)
        spectrogram_seq_len = size(X_test_snr{j}, 2);
        original_hop_sequence = Y_test_snr_raw{j};
        aligned_labels = zeros(1, spectrogram_seq_len);
        for t = 1:spectrogram_seq_len
            center_sample_index = (t - 1) * stft_hop_size + floor(stft_params.windowLength / 2);
            hop_index_for_this_frame = floor((center_sample_index - 1) / samples_per_hop) + 1;
            hop_index_for_this_frame = max(1, min(hop_index_for_this_frame, simParams.numHops));
            aligned_labels(t) = original_hop_sequence(hop_index_for_this_frame);
        end
        Y_test_aligned{j} = aligned_labels;
    end
    
    numClasses = 2^3; % k=3, so 8 classes (0-7)
    classNames = string(0:(numClasses-1));
    Y_test_cat = cell(size(Y_test_aligned));
    for j = 1:length(Y_test_aligned)
        Y_test_cat{j} = categorical(Y_test_aligned{j}, 0:(numClasses-1), classNames);
    end

    % c. Classify using the network
    fprintf('Classifying %d sequences...\n', length(X_test_snr));
    Y_pred = classify(net, X_test_snr, 'MiniBatchSize', 32);

    % --- NEW: Display example predictions ---
    % Let's display a few examples for a specific SNR, e.g., 0 dB
    if snr_db == 0
        fprintf('--- Example Predictions for SNR = 0 dB ---\n');
        num_examples_to_show = 3; % Reduced for better readability
        num_hops_to_display = 20; % Reduced for better readability
        for ex_idx = 1:min(num_examples_to_show, length(Y_pred))
            fprintf('Example %d:\n', ex_idx);

            % Convert categorical back to numeric index (0-based)
            gt_indices = double(string(Y_test_cat{ex_idx}(1:num_hops_to_display)));
            pred_indices = double(string(Y_pred{ex_idx}(1:num_hops_to_display)));

            % Map indices to frequencies (MHz for readability)
            gt_freqs_mhz = hopset(gt_indices + 1) / 1e6;
            pred_freqs_mhz = hopset(pred_indices + 1) / 1e6;

            % Create formatted strings for display
            ground_truth_str = sprintf('%.2f ', gt_freqs_mhz);
            predicted_str = sprintf('%.2f ', pred_freqs_mhz);

            fprintf('  Ground Truth (MHz): %s...\n', ground_truth_str);
            fprintf('  Predicted (MHz):    %s...\n', predicted_str);
            
            % Highlight differences
            is_different = Y_test_cat{ex_idx}(1:num_hops_to_display) ~= Y_pred{ex_idx}(1:num_hops_to_display);
            diff_indices = find(is_different);
            if isempty(diff_indices)
                fprintf('  (No errors in the first %d hops)\n\n', num_hops_to_display);
            else
                fprintf('  (Errors at hop indices: %s)\n\n', num2str(diff_indices'));
            end
        end
    end
    % --- END NEW SECTION ---

    % d. Calculate accuracy
    num_correct = 0;
    num_total = 0;
    for j = 1:length(Y_pred)
        % Compare predicted sequence with the ground truth sequence
        num_correct = num_correct + sum(Y_pred{j} == Y_test_cat{j});
        num_total = num_total + numel(Y_pred{j});
    end
    accuracy = num_correct / num_total;
    accuracies(i) = accuracy;
    
    fprintf('Accuracy @ %d dB = %.4f\n', snr_db, accuracy);
end

% --- 5. Plot Results ---
fprintf('Plotting results...\n');
figure;
plot(snr_levels_db, accuracies * 100, 'o-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Hop Prediction Accuracy (%)');
title('LSTM Model Performance vs. SNR');
ylim([0 100]);
set(gca, 'FontSize', 12);

fprintf('Evaluation complete.\n');
