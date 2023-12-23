close all 
clear all

%% Load final results for figure

% Large scouts PLDC:
load('final_mats_for_figures/large_scouts_single_run.mat')
large_scouts_params

% % Small scouts PLDC:
% load('final_mats_for_figures/small_scouts_single_run.mat')
% small_scouts_params

%% Plot spatial delay pattern

line_width = 3;

hCbarYlabel = 'Delay [ms]';
hCbarYlablFontSize = 14;

signal_color = ['g','k'];

%load channel positions
path_to_channel_positions = [get_analysis_code_base_path() 'spatial_data/channel_locations.mat'];
electrode_channel_poitions = load(path_to_channel_positions );
electrode_channel_poitions = electrode_channel_poitions.channel_locations;

if isempty(plotting_config.start_end_wave)
    analysis_start_end_wave = [1 size(single_run.electrode_set.electrode_signals,2)];
end
[phase_event_heights, phase_event_locs] = single_run.get_phase_events("first_max",analysis_start_end_wave );
t = [1:size(single_run.electrode_set.electrode_signals,2)] / plotting_config.fs;

f1 = figure("OuterPosition",[680,558,600,420]);
ax = axes("Position",[0.13,0.11,0.7,0.815]);
ax.PositionConstraint = "innerposition";
single_run.signal_source.plot_source_positions(gca);
hold on
if plot_flipped_scouts
    flip_around = mean(electrode_channel_poitions(single_run.electrode_grid(3,:),2));
    single_run.signal_source.plot_source_positions(gca,true,flip_around);
end
scatter3(electrode_channel_poitions(single_run.electrode_subset,1),electrode_channel_poitions(single_run.electrode_subset,2),electrode_channel_poitions(single_run.electrode_subset,3),600,1000 * (phase_event_locs - min(phase_event_locs)) / plotting_config.fs,'filled','MarkerEdgeColor',[0 0 0])
axis off
hCbar = colorbar;
c_bar_ticks = get(hCbar,'YTick');
set(hCbar,{'YTick','YTickLabel'},{[c_bar_ticks(1),c_bar_ticks(end)],[c_bar_ticks(1),c_bar_ticks(end)]})
set(hCbar,'Position',[0.85967177743371,0.769757414546963,0.016903307680488,0.144549763033177])
ylabel(hCbar,hCbarYlabel,{'FontSize','Position'},{8,[3.506666618982955,7.900007534027101,0]});


%% Choose relevant params
% large_scouts_params
small_scouts_params
% large_scouts_filtered_params

%% Load and calculate spatial pattern

%Channel Positions
path_to_channel_positions = [get_analysis_code_base_path() 'spatial_data/channel_locations.mat'];

% Apply path configurations and load
[electrode_channel_poitions, cortex, simulation_table, electrode_grid, save_path] = apply_path_configurations_and_load(path_to_channel_positions,cortex_name,path_to_results_table,electrode_grid_name, ''); 
simulation_table = simulation_table.two_scouts_results;

%%% Get relevant row and plot head
row_num = get_row_num_by_dt_dx(simulation_table,dt,dx);
row = simulation_table(row_num ,:);

%%% Create cortex, electrode and scout objects
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
electrode_set = ElectrodeSet(electrode_signals(:,plotting_config.crop_at(1):plotting_config.crop_at(2)), electrode_channel_poitions,plotting_config);
sources = Scouts.init_source_from_table_row(row, cortex);
single_run = SimulationAnalysis(sources, electrode_set, electrode_grid, plotting_config.fs, cortex, dt, dx);

single_run.electrode_set.parent_analysis = single_run;
single_run.signal_source.parent_analysis = single_run;

% save('final_mats_for_figures/single_run.mat',"single_run")
