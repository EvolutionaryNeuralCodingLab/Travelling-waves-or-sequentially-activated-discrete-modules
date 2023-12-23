clear all
close all

%% Plot filtered vs non-filtered electorde signal

load('final_mats_for_figures/filtered_vs_non_filtered_data.mat','t','channels_for_signal','unfiltered_signal','filtered_signal')

line_width = 3;

f1 = figure;
plot(t,unfiltered_signal,'LineWidth',line_width)
hold on
plot(t,filtered_signal,'LineWidth',line_width)
legend('Unfiltered signal','Filtered signal')
xlabel('Time [s]')
ylabel('Electrode Simulated Signal [AU]')

%% Calculate filtered vs non-filtered electorde signal

large_scouts_filtered_params

%Channel Positions
path_to_channel_positions = [get_analysis_code_base_path() 'spatial_data/channel_locations.mat'];

% Apply path configurations and load
[electrode_channel_poitions, cortex, simulation_table, electrode_grid, save_path] = apply_path_configurations_and_load(path_to_channel_positions,cortex_name,path_to_results_table,electrode_grid_name, ''); 
simulation_table = simulation_table.two_scouts_results;

% Get relevant row and plot head
row_num = get_row_num_by_dt_dx(simulation_table,dt,dx);
row = simulation_table(row_num ,:);

if plotting_config.duplicate_hemisphere
   temp_signal = row.EEG_recordings{1};
   flipped_electrode_grid = flipud(electrode_grid);
   temp_signal_flipped = temp_signal;
   temp_signal_flipped(electrode_grid(:),:) = temp_signal_flipped(flipped_electrode_grid(:),:);
   temp_signal_duplicated = temp_signal;
   temp_signal_duplicated(electrode_grid(:),:) = temp_signal(electrode_grid(:),:) + temp_signal_flipped(electrode_grid(:),:);
   electrode_signals = temp_signal_duplicated;
else
   electrode_signals = row.EEG_recordings{1};
end

%%% Unfiltered signal
plotting_config.ws = [];
plotting_config = fill_missing_configs_with_defaults(plotting_config);
electrode_set_unfiltered = ElectrodeSet(electrode_signals(:,plotting_config.crop_at(1):plotting_config.crop_at(2)), electrode_channel_poitions,plotting_config);

%%% Filtered signal
plotting_config.ws = [10 12];
plotting_config.flip_signal = 0;
plotting_config.pad_samples = 1000;
plotting_config.n_reps = 100;
electrode_set_filtered = ElectrodeSet(electrode_signals(:,plotting_config.crop_at(1):plotting_config.crop_at(2)), electrode_channel_poitions,plotting_config);

channels_for_signal = 53;
t = (1:size(electrode_set_unfiltered.electrode_signals,2)) / plotting_config.fs;
unfiltered_signal = electrode_set_unfiltered.electrode_signals(channels_for_signal,:);
filtered_signal = electrode_set_filtered.electrode_signals(channels_for_signal,:);

save('final_mats_for_figures/filtered_vs_non_filtered_data.mat','t','channels_for_signal','unfiltered_signal','filtered_signal')
