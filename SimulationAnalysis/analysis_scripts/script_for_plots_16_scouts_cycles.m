%% Plot elecrode map with signal

clear all
close all

%% Definitions and paths

dt = 0.2;
dx = 6;
% plotting_config.ws = 1 ./ [1.6 1.4];
plotting_config.ws = [10 12];
% plotting_config.fs = 600;
plotting_config.fs = 10000;
plotting_config.start_end_wave = []; % [500, 1500];
plotting_config.crop_at = [260, 5000];
plotting_config.pad_samples = 1000;
plotting_config.flip_signal = 0;
plotting_config.n_reps = 100;
plotting_config.duplicate_hemisphere = 1;
plotting_config = fill_missing_configs_with_defaults(plotting_config);

electrode_grid_name = 'En7_lr_symmetric';
cortex_name = "306716V";

results_filename = '/media/sil2/Data/Yuval O/Rotem Simulations/20230702_RightHemisphereSmallScouts_306716/simulation results.mat';

path_to_channel_positions = '/media/sil2/Data/Yuval O/Rotem Simulations/SimulationAnalysis/spatial_data/channel_locations.mat';
figs_dir = '16_scouts_cycles/';

wave_metrics_matrix_path = '/media/sil2/Data/Yuval O/Rotem Simulations/20230702_RightHemisphereSmallScouts_306716/figures/flip_signal/En7_lr_symmetric/Two Scouts - wave metrics - first maxima.mat';



%% load stuff

wave_metrics_matrix = load(wave_metrics_matrix_path);
%electrode positions
electrode_channel_poitions = load(path_to_channel_positions);
electrode_channel_poitions = electrode_channel_poitions.channel_locations;
%cortex vertices
cortex = load_cortex(cortex_name);
%simulation table
simulation_table = load(results_filename);
simulation_table = simulation_table.two_scouts_results;
%electrode grid
electrode_grids_path = '/media/sil2/Data/Yuval O/Rotem Simulations/SimulationAnalysis/spatial_data/elecrode_grids/'; %TODO: move to consts
electrode_grid_path = [electrode_grids_path,electrode_grid_name,'.txt'];
electrode_grid = readmatrix(electrode_grid_path);

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

disp(['Selected simulation''s PLDC: ' num2str(wave_metrics_matrix.PLDCs(find(wave_metrics_matrix.delta_Ts==dt),find(wave_metrics_matrix.delta_Xs==dx)))])

%% Plot Electrodes' time to first max - 16 scouts cycle
close all
save_path = '/media/sil2/Data/Yuval O/Rotem Simulations/figures_for_paper/16_scouts_cycles/En7_lr_symmetric/';
line_width = 3;
channels_for_signal = [305,283];

% Plot PLDC
f = figure("OuterPosition",[680,558,580,420]);
ax = axes("Position",[0.13,0.11,0.6,0.815]);
ax.PositionConstraint = "innerposition";
T_start_end = [1,length(wave_metrics_matrix.delta_Ts)];
X_start_end = [2,length(wave_metrics_matrix.delta_Xs)];
img_range = [0.85,1];
f = plot_phase_space(wave_metrics_matrix.PLDCs,wave_metrics_matrix.delta_Ts, wave_metrics_matrix.delta_Xs, T_start_end,X_start_end,0,img_range,ax);
hCbar = colorbar;
c_bar_ticks = get(hCbar,'YTick');
set(hCbar,{'YTick','YTickLabel'},{[c_bar_ticks(1),c_bar_ticks(end)],[c_bar_ticks(1),c_bar_ticks(end)]})
set(hCbar,'Position',[0.745652728932942,0.778659491698299,0.016903307680488,0.144549763033177])
ylabel(hCbar,'Correlation',{'FontSize','Position'},{8,[4.606666586200397,0.925000071525574,0]});
set(gcf,'PaperPositionMode','auto');%print name -dpdf -painters;
saveJpegAndFig(f,save_path,'PLDC',1,'saveEPS',1)


%Plot Spatial delay pattern

start_end_wave = [];

hCbarYlabel = 'Delay [ms]';
hCbarYlablFontSize = 14;

if isempty(plotting_config.start_end_wave)
    analysis_start_end_wave = [1 size(single_run.electrode_set.electrode_signals,2)];
end
[phase_event_heights, phase_event_locs] = single_run.get_phase_events("first_max",analysis_start_end_wave );
t = (1:size(electrode_set.electrode_signals,2)) / plotting_config.fs;

f1 = figure("OuterPosition",[680,558,600,420]);
ax = axes("Position",[0.13,0.11,0.7,0.815]);
ax.PositionConstraint = "innerposition";sources.plot_source_positions(gca);
hold on
scatter3(electrode_channel_poitions(single_run.electrode_subset,1),electrode_channel_poitions(single_run.electrode_subset,2),electrode_channel_poitions(single_run.electrode_subset,3),600,1000 * (phase_event_locs - min(phase_event_locs)) / plotting_config.fs,'filled','MarkerEdgeColor',[0 0 0])
scatter3(electrode_channel_poitions(channels_for_signal(1),1),electrode_channel_poitions(channels_for_signal(1),2),electrode_channel_poitions(channels_for_signal(1),3),[],signal_color(1), ...
    'filled')
scatter3(electrode_channel_poitions(channels_for_signal(2),1),electrode_channel_poitions(channels_for_signal(2),2),electrode_channel_poitions(channels_for_signal(2),3),[],signal_color(2),'filled')
axis off
hCbar = colorbar;
c_bar_ticks = get(hCbar,'YTick');
set(hCbar,{'YTick','YTickLabel'},{[c_bar_ticks(1),c_bar_ticks(end)],[c_bar_ticks(1),c_bar_ticks(end)]})
set(hCbar,'Position',[0.904822279105951,0.760855337395628,0.016903307680488,0.144549763033177])
ylabel(hCbar,hCbarYlabel,{'FontSize','Position'},{8,[3.506666618982955,6.450006151199341,0]});
saveJpegAndFig(f1,save_path,'spatial_lag',1,'saveEPS',1)


f2 = figure;
t = (1:size(electrode_set.electrode_signals,2)) / plotting_config.fs;
signal_color = ['g','k'];
plot(t, 10e6*(electrode_set.electrode_signals(channels_for_signal(1),:) + max(electrode_set.electrode_signals(channels_for_signal(1),:))),signal_color(1),'LineWidth',line_width)
hold on
plot(t, 10e6*electrode_set.electrode_signals(channels_for_signal(2),:),signal_color(2),'LineWidth',line_width)
xlabel('Time [s]')
ylabel('Electrode Simulated Signal [AU]')
set(gcf,'PaperPositionMode','auto');%print name -dpdf -painters;
saveJpegAndFig(f2,save_path,'electrode_signal',1,'saveEPS',1)

f3 = figure;
t = (1:size(single_run.signal_source.signal(1,:),2)) / plotting_config.fs;
plot(t, single_run.signal_source.signal(1,:) + max(single_run.signal_source.signal(2,:)),'b','LineWidth',line_width)
hold on
plot(t, single_run.signal_source.signal(2,:),'r','LineWidth',line_width)
xlabel('Time [s]')
ylabel('Scouts Signal [AU]')
xlim([0 0.5])
set(gcf,'PaperPositionMode','auto');%print name -dpdf -painters;
saveJpegAndFig(f3,save_path,'scout_signal',1,'saveEPS',1)

% set(gca,'YTick',[])

%% Plot Electrodes' time to first max - 16 scouts cycle
close all
save_path = '/media/sil2/Data/Yuval O/Rotem Simulations/figures_for_paper/16_scouts_cycles/En7/';

% Plot PLDC
f = figure("OuterPosition",[680,558,580,420]);
ax = axes("Position",[0.13,0.11,0.6,0.815]);
ax.PositionConstraint = "innerposition";
T_start_end = [1,length(wave_metrics_matrix.delta_Ts)];
X_start_end = [2,length(wave_metrics_matrix.delta_Xs)];
img_range = [0.85,1];
f = plot_phase_space(PLDCs,delta_Ts, delta_Xs, T_start_end,X_start_end,0,img_range,ax);
hCbar = colorbar;
c_bar_ticks = get(hCbar,'YTick');
set(hCbar,{'YTick','YTickLabel'},{[c_bar_ticks(1),c_bar_ticks(end)],[c_bar_ticks(1),c_bar_ticks(end)]})
set(hCbar,'Position',[0.745652728932942,0.778659491698299,0.016903307680488,0.144549763033177])
ylabel(hCbar,'Correlation',{'FontSize','Position'},{8,[4.606666586200397,0.925000071525574,0]});
set(gcf,'PaperPositionMode','auto');%print name -dpdf -painters;
saveJpegAndFig(f,save_path,'PLDC',1,'saveEPS',1)


%Plot Spatial delay pattern
start_end_wave = [];
channels_for_signal = [305,283];
line_width = 3;

hCbarYlabel = 'Delay [ms]';
hCbarYlablFontSize = 14;

signal_color = ['g','k'];
if isempty(plotting_config.start_end_wave)
    analysis_start_end_wave = [1 size(single_run.electrode_set.electrode_signals,2)];
end
[phase_event_heights, phase_event_locs] = single_run.get_phase_events("first_max",analysis_start_end_wave );
t = (1:size(electrode_set.electrode_signals,2)) / plotting_config.fs;

f1 = figure;
sources.plot_source_positions(gca);
hold on
flip_around = mean(electrode_channel_poitions(electrode_grid(3,:),2));
sources.plot_source_positions(gca,true,flip_around);
scatter3(electrode_channel_poitions(single_run.electrode_subset,1),electrode_channel_poitions(single_run.electrode_subset,2),electrode_channel_poitions(single_run.electrode_subset,3),600,1000 * (phase_event_locs - min(phase_event_locs)) / plotting_config.fs,'filled','MarkerEdgeColor',[0 0 0])
scatter3(electrode_channel_poitions(channels_for_signal(1),1),electrode_channel_poitions(channels_for_signal(1),2),electrode_channel_poitions(channels_for_signal(1),3),[],signal_color(1), ...
    'filled')
scatter3(electrode_channel_poitions(channels_for_signal(2),1),electrode_channel_poitions(channels_for_signal(2),2),electrode_channel_poitions(channels_for_signal(2),3),[],signal_color(2),'filled')
axis off
hCbar = colorbar;
c_bar_ticks = get(hCbar,'YTick');
set(hCbar,{'YTick','YTickLabel'},{[c_bar_ticks(1),c_bar_ticks(end)],[c_bar_ticks(1),c_bar_ticks(end)]})
set(hCbar,'Position',[0.904822279105951,0.760855337395628,0.016903307680488,0.144549763033177])
ylabel(hCbar,hCbarYlabel,{'FontSize','Position'},{10,[3.995959546830919,108.0883361399174,0]});
cd('/media/sil2/Data/Yuval O/Rotem Simulations/figures_for_paper/16_scouts_cycles/En7_lr_symmetric/')
% set(gcf,'PaperPositionMode','auto');print spatial_lag_ghost_scouts_gray -dpdf -painters;
% set(gcf,'PaperPositionMode','auto');print spatial_lag_ghost_scouts_gray_filled -dpdf -painters;
set(gcf,'PaperPositionMode','auto');print spatial_lag_ghost_scouts_colored_filled -dpdf -painters;
% set(gcf,'PaperPositionMode','auto');print spatial_lag_ghost_scouts_colored_hollowed -dpdf -painters;
close(f1)

f2 = figure;
plot(t, 10e6*(electrode_set.electrode_signals(channels_for_signal(1),:) + max(electrode_set.electrode_signals(channels_for_signal(1),:))),signal_color(1),'LineWidth',line_width)
hold on
plot(t, 10e6*electrode_set.electrode_signals(channels_for_signal(2),:),signal_color(2),'LineWidth',line_width)
xlabel('Time [s]')
ylabel('Electrode Simulated Signal [AU]')

f3 = figure;
t = (1:size(single_run.signal_source.signal(1,:),2)) / plotting_config.fs;
plot(t, single_run.signal_source.signal(1,:) + max(single_run.signal_source.signal(2,:)),'b','LineWidth',line_width)
hold on
plot(t, single_run.signal_source.signal(2,:),'r','LineWidth',line_width)
xlabel('Time [s]')
ylabel('Scouts Signal [AU]')

%% OLD BACKUP SCRIPT - 
% % Plot Electrodes' signal - 16 scouts cycles
% 
% close all
% 
% half_point = 725;
% f = figure;
% sources.plot_source_positions(gca);
% colorbar
% hold on
% scatter3(electrode_channel_poitions(:,1),electrode_channel_poitions(:,2),electrode_channel_poitions(:,3),200,electrode_set.electrode_signals(:,half_point)*10e7,'filled')

axis off