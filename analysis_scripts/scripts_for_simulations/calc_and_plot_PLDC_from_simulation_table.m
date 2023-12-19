clear all
close all


%% Uncomment relevant params

large_scouts_params % Main simulation used
% small_scouts_params
% large_scouts_filtered_params

%% Load parameters (paths inside)

%Channel Positions
path_to_channel_positions = [get_analysis_code_base_path() 'spatial_data/channel_locations.mat'];

%Base save path (subset grid subfolder will be added automatically
save_path = ''; % Leaing empty prevents saving metrics; if any plotting option is true it'll save it in current dir

% Apply path configurations and load
[electrode_channel_poitions, cortex, simulation_table, electrode_grid, save_path] = apply_path_configurations_and_load(path_to_channel_positions,cortex_name,path_to_results_table,electrode_grid_name, save_path); 
simulation_table = simulation_table.two_scouts_results;


%% Load results and run analysis
rng(0) % Set seed for reproducibility
[dip_p_values, PLDC, phase_space_dT, phase_space_dX, delta_Ts, delta_Xs] = analyze_simulations_from_table(simulation_table, electrode_channel_poitions, electrode_grid, source_type, sampling_rate, cortex, plotting_config, save_path, files_prefix);

%% Plot PLDC

[all_scouts_dipoles_positions, mean_electrode_dist, scouts_dist_per_electrode_dist, scouts_stds, scouts_dist_per_std] = ...
    calc_distances(dxs_for_distance, dt, electrode_grid, cortex, simulation_table, electrode_channel_poitions);

%PLDC vs dT - legend is dx
f=figure;
plot(delta_Ts,PLDC)
ylim(plotting_config.PLDC_summary_scale)
ylabel('PLDC')
xlabel('\DeltaT / \sigma_T')
legend(num2str(round(scouts_dist_per_std,1)))
