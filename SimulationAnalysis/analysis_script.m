%Use this script to set relevant paths (cortex,simulation results, save
%path, etc) and run relevant analysis function (e.g.
%analyze_simulations_from_table)

clear all
close all

%% Set paths and var names
%Channel Positions
path_to_channel_positions = '/media/sil2/Data/Yuval O/Rotem Simulations/SimulationAnalysis/spatial_data/channel_locations.mat';

%Cortex Name
% cortex_name = "306716V";
cortex_name = "15002V";

%Simulation table
results_dir = '/media/sil2/Data/Yuval O/Rotem Simulations/six_large_scouts/';
results_filename = 'simulation results.mat';
path_to_results_table = [results_dir, results_filename];

%Base save path (subset grid subfolder will be added automatically
save_path = [results_dir 'figures/'];

%Choose electrode subset grid
electrode_grid_name = 'En';

% Apply path configurations and load
[electrode_channel_poitions, cortex, simulation_table, electrode_grid, save_path] = apply_path_configurations_and_load(path_to_channel_positions,cortex_name,path_to_results_table,electrode_grid_name, save_path); 

%set sampling rate
sampling_rate = 600; %Hz

%%TODO: Add plot to see how En and sources look like toghat

%% Set plotting options

plotting_config.with_plots = 1; %Plot and save all possible plots
plotting_config.plot_summary = 0; %Just plot the final dip and PLDC phase spaces
plotting_config.plot_signal = 0; %Plots and save recorded signals from electrode
plotting_config.plot_source_signal = 0; %Plot and save source signals
plotting_config.plot_source_position = 0; %Plot and save source positions
plotting_config.duplicate_hemisphere = 0; %flip the signal o f electrode grid up and add them to original signal (used for when onyl one hemisphere is activated)
% plotting_config.percentile = 0.05; %When analyzing results, drop top and bottom plotting_config.percentile of crossing times
plotting_config.percentile = 0;
plotting_config.interpolate = 0; %Add additional crossings between any two crossings for smoothness
plotting_config.ignore_source = 0;
% plotting_config.ws = pi ./ [2200 800];
plotting_config.ws = 1 ./ [1.6 1.4];
plotting_config.fs = sampling_rate;
plotting_config.start_end_wave = [500 1500];

%% Run this for scouts scripts

simulation_table = simulation_table.two_scouts_results;
source_type = 'scouts';
files_prefix = "Two Scouts - ";

%% Run this for dipole script

%% Start analysis

[dip_p_values, PLDC, phase_space_dT, phase_space_dX] = analyze_simulations_from_table(simulation_table, electrode_channel_poitions, electrode_grid, source_type, sampling_rate, cortex, plotting_config, save_path, files_prefix);
