clear all
close all

%% Load parameters (paths inside)

% Choose relevant params
small_scouts_no_filtering_params
% small_scouts_with_filtering_params

%Channel Positions
path_to_channel_positions = [get_analysis_code_base_path() 'spatial_data/channel_locations.mat'];

%Base save path (subset grid subfolder will be added automatically
save_path = [results_dir 'figures/'];

% Apply path configurations and load
[electrode_channel_poitions, cortex, simulation_table, electrode_grid, save_path] = apply_path_configurations_and_load(path_to_channel_positions,cortex_name,path_to_results_table,electrode_grid_name, save_path); 


%% Load results and run analysis

simulation_table = simulation_table.two_scouts_results;
[dip_p_values, PLDC, phase_space_dT, phase_space_dX] = analyze_simulations_from_table(simulation_table, electrode_channel_poitions, electrode_grid, source_type, sampling_rate, cortex, plotting_config, save_path, files_prefix);
