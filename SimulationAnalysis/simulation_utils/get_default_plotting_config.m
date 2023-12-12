function default_plotting_config = get_default_plotting_config()
%GET_DEFAULT_default_plotting_config Summary of this function goes here
%   Detailed explanation goes here

default_plotting_config.with_plots = 1; %Plot and save all possible plots
default_plotting_config.plot_summary = 0; %Just plot the final dip and PLDC phase spaces
default_plotting_config.plot_signal = 0; %Plots and save recorded signals from electrode
default_plotting_config.plot_source_signal = 0; %Plot and save source signals
default_plotting_config.plot_source_position = 0; %Plot and save source positions
default_plotting_config.duplicate_hemisphere = 0; %flip the signal of electrode grid up and add them to original signal (used for when onyl one hemisphere is activated)
default_plotting_config.percentile = 0; %When analyzing results, drop top and bottom default_plotting_config.percentile of crossing times
default_plotting_config.interpolate = 0; %Add additional crossings between any two crossings for smoothness
default_plotting_config.ignore_source = 0;
default_plotting_config.start_end_wave = [];% [500, 1500]; %samples in between which to start looking at phase events
default_plotting_config.fs = 600; %sampling rate
default_plotting_config.noise_gauss_params = []; %[gausian_mean, gausian_std]
default_plotting_config.save_analysis_objects = 0;
% bandpassing parameters
%%%Before bandpassing, pad_samples will be repeated n_reps times starting
%%%from crop_at(1). Same for the end
default_plotting_config.ws = [];% [0.5 2];%1 ./ [1.6 1.4]; %bandpass frequencies
default_plotting_config.crop_at = [];% [260, 5000]; 
default_plotting_config.pad_samples = 0; %1000;
default_plotting_config.flip_signal = 0;
default_plotting_config.n_reps = 100;

