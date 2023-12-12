% Path to simulation table
results_dir = '/media/sil2/Data/Yuval O/Rotem Simulations/20230702_RightHemisphereSmallScouts_306716/';
results_filename = 'simulation results.mat';
path_to_results_table = [results_dir, results_filename];

cortex_name = "306716V";
electrode_grid_name = 'electrode_grid_for_small_scouts';
sampling_rate = 600; %Hz 
source_type = 'scouts';
files_prefix = "Two Scouts - "; % Saved files prefix

plotting_config.duplicate_hemisphere = 1; %flip the signal of electrode grid up and add them to original signal (used for when onyl one hemisphere is activated)
plotting_config.save_analysis_objects = 0;

plotting_config.percentile = 0; % When analyzing results, drop top and bottom plotting_config.percentile of crossing times
plotting_config.interpolate = 0; % Add additional crossings between any two crossings for smoothness
plotting_config.ignore_source = 0; % Some plots rely on having a defined signal source. If a run is done without one - use this to avoid bugs
plotting_config.fs = 10000; % Sampling rate (Hz)
plotting_config.start_end_wave = []; % Look at phase crossings only within these samples
plotting_config.crop_at = [260, 5000];
plotting_config.pad_samples = 1000;
plotting_config.n_reps = 100;
plotting_config.flip_signal = 0;

plotting_config.with_plots = 0; % Plot and save all possible plots. 
% For the following plots, if with_plots is false, you can set specific
% plots. If with_plots true, they will be plotted even if 0.
plotting_config.plot_summary = 0; % Just plot the final dip and PLDC phase spaces
plotting_config.plot_signal = 0; % Plots and save recorded signals from electrode
plotting_config.PLDC_summary_scale = [0.85 1];
plotting_config.plot_source_signal = 0; % Plot and save source signals
plotting_config.plot_source_position = 0; %P lot and save source positions