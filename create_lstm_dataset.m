% --- create_lstm_dataset.m ---
% 
% Generates a dataset of FHSS signals with AWGN for use in an LSTM model.
% The dataset will contain signals with varying hopsets and SNR levels.

% --- Configuration ---

% SNR levels in dB
snr_levels_db = -10:2:10;

% Number of different hopsets to generate
num_hopsets = 10;

% Number of signals to generate per hopset and SNR level
num_signals_per_config = 100;

% --- Dataset Initialization ---
X_data = {}; % Cell array for input sequences (time-frequency data)
Y_data = {}; % Cell array for output labels (hop sequences)

% --- Data Generation Loop ---
fprintf('Starting dataset generation...\n');

for i = 1:length(snr_levels_db)
    snr_db = snr_levels_db(i);
    fprintf('Processing SNR = %d dB\n', snr_db);

    for j = 1:num_hopsets
        fprintf('  Generating data for hopset %d/%d\n', j, num_hopsets);
        
        % In a real implementation, we would modify the hop generation here.
        % For now, we'll assume the default hop generation in the transmitter
        % is sufficient and we'll just generate multiple realizations.
        % A more advanced version would change the Gold code polynomials.
        
        for k = 1:num_signals_per_config
            % 1. Run the transceiver simulation with AWGN
            % This part needs to be adapted from main_02_run_transceiver_awgn.m
            % We will call a modified simulation function.
            
            % [received_signal, original_hops] = run_fhss_simulation(snr_db, hopset_seed);
            
            % 2. Perform Time-Frequency Analysis to get LSTM input
            % This will use functions from src/analysis/
            
            % tfr_data = compute_tfr(received_signal);
            
            % 3. Store data
            % X_data{end+1} = tfr_data;
            % Y_data{end+1} = original_hops;
            
        end
    end
end

% --- Save Dataset ---
fprintf('Saving dataset to file...\n');
% save('data/synthetic/lstm_dataset_awgn.mat', 'X_data', 'Y_data', '-v7.3');

fprintf('Dataset generation complete.\n');
