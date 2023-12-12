clear all
close all

%% Load parameters (paths inside)

% Choose relevant params
% large_scouts_params
large_scouts_filtered_params
% small_scouts_params

%Channel Positions
path_to_channel_positions = [get_analysis_code_base_path() 'spatial_data/channel_locations.mat'];

%Base save path (subset grid subfolder will be added automatically
save_path = ['/media/sil2/Data/Yuval O/Rotem Simulations/New repo tests/figures/filtered/'];

% Apply path configurations and load
[electrode_channel_poitions, cortex, simulation_table, electrode_grid, save_path] = apply_path_configurations_and_load(path_to_channel_positions,cortex_name,path_to_results_table,electrode_grid_name, save_path); 


%% Load results and run analysis
rng(0) % Set seed for reproducibility
simulation_table = simulation_table.two_scouts_results;
[dip_p_values, PLDC, phase_space_dT, phase_space_dX, delta_Ts, delta_Xs] = analyze_simulations_from_table(simulation_table, electrode_channel_poitions, electrode_grid, source_type, sampling_rate, cortex, plotting_config, save_path, files_prefix);

%% Plot PLDC

% wave_metrics_matrix_path = '/media/sil2/Data/Yuval O/Rotem Simulations/221204 simulations/large size scouts/multiple bumps/figures/En9/Two Scouts - wave metrics - first maxima.mat';
% wave_metrics_matrix = load(wave_metrics_matrix_path);
%distances calculated in /media/sil2/Data/Yuval O/Rotem Simulations/figures_for_paper/spatial_calcs.m
distances = load('/media/sil2/Data/Yuval O/Rotem Simulations/figures_for_paper/spatial_calcs/scotus_electrode_distances_large_scouts.mat');
% save_path = '/media/sil2/Data/Yuval O/Rotem Simulations/figures_for_paper/PLDC_as_lines/';

%PLDC vs dT - legend is dx
f=figure;
plot(wave_metrics_matrix.delta_Ts,wave_metrics_matrix.PLDCs)
% legend(num2str(round(distances.scouts_mean_distances*100,2,'significant')))
legend(num2str(round(distances.scouts_mean_distances*100,1)))
ylim(plotting_config.PLDC_summary_scale)
ylabel('PLDC')
xlabel('\DeltaT / \sigma_T')
legend(num2str(round(distances.scouts_dist_per_std,1)))
