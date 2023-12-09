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

%%TODO: Add plot to see how En and sources look like toghat

%% Set plotting options

plotting_config.with_plots = 1; %Plot and save all possible plots
plotting_config.plot_summary = 0; %Just plot the final dip and PLDC phase spaces
plotting_config.plot_signal = 0; %Plots and save recorded signals from electrode
plotting_config.plot_source_signal = 0; %Plot and save source signals
plotting_config.plot_source_position = 0; %Plot and save source positions
plotting_config.duplicate_hemisphere = 0; %flip the signal of electrode grid up and add them to original signal (used for when onyl one hemisphere is activated)
% plotting_config.percentile = 0.05; %When analyzing results, drop top and bottom plotting_config.percentile of crossing times
plotting_config.percentile = 0;
plotting_config.interpolate = 0; %Add additional crossings between any two crossings for smoothness

%% Run this for scouts scripts

simulation_table = simulation_table.two_scouts_results;
source_type = 'scouts';
files_prefix = "Two Scouts - ";

%% Run this for dipole script

%% Start analysis

[dip_p_values, PLDC, phase_space_dT, phase_space_dX] = analyze_simulations_from_table(simulation_table, electrode_channel_poitions, electrode_grid, source_type, cortex, plotting_config, save_path, files_prefix);

%%



%extract small electrode grid to focus analysis on
% electrode_map = readmatrix('/media/sil2/Data/Yuval O/Rotem Simulations/spatial_skull_info/En.txt');
electrode_map = readmatrix('/media/sil2/Data/Yuval O/Rotem Simulations/spatial_skull_info/En_smaller3.txt');
% electrode_map = readmatrix('/media/sil2/Data/Yuval O/Rotem Simulations/EEG2scouts_grid 2/EEG2scouts_grid/En_smaller2.txt');
% electrode_map = readmatrix('/media/sil2/Data/Yuval O/Rotem Simulations/EEG2scouts_grid 2/EEG2scouts_grid/En2.txt');
% electrode_map = importdata('/media/sil2/Data/Yuval O/Rotem Simulations/multiple scouts/En_two_row_opt2.txt');


%%
% Yuval, what are these?
delta_Ts = unique(T.deltaT);
delta_Xs = unique(T.distance);

%%




n_rows = height(T);
n_scouts = 2;

cors = zeros(length(delta_Ts),length(delta_Xs));
dip_p_values = cors;

elecrode_nums = electrode_map(~isnan(electrode_map));
nCh_sub = length(elecrode_nums);

electrode_columns = repmat(1:size(electrode_map,2),size(electrode_map,1),1);
n_electrode_columns = size(electrode_map,2);
column_colors = parula(n_electrode_columns);

% channel coordinates
x = channel_locations(:,1);
y = channel_locations(:,2);
z = channel_locations(:,3);        
x_sub = x(elecrode_nums);
y_sub = y(elecrode_nums);
z_sub = z(elecrode_nums);   
 
t_sanity = cors;
x_sanity = cors;
    
for i=1:n_rows
   row = T(i,:);
   delta_T = row.deltaT;
   delta_X = row.distance;
   
   simulated_recording = row.EEG_recordings{1,1};   
 
   %plot 3d coordinates all electrodes
    if (with_plots || plot_scout_position) && delta_T == delta_Ts(1)
        f=figure;
        text(x,y,z, string(1:numel(x)));
        hold on
        scout_colors = 'rb';
        for scout_num=1:n_scouts
            scout_vertices_inds = row.scouts{1,1}{1,scout_num}.Vertices;
%             scout_vertices_locs = cortex_vars.cortex306716V.Vertices(scout_vertices_inds,:);
            scout_vertices_locs = cortex_vars.cortex.Vertices(scout_vertices_inds,:);
            scatter3(scout_vertices_locs(:,1),scout_vertices_locs(:,2),scout_vertices_locs(:,3),scout_colors(scout_num),'filled');
        end
        % There are two channels situated at the same location: 60,337

        axis([min(x), max(x), min(y), max(y),min(z), max(z)])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('Elecrode Position (scouts vertices in color) - 3d')
        saveJpegAndFig(f,save_path,char(files_prefix + "- electode position 3d - Delta_X " + delta_X),1)
        close(f)
    end
    
    if (with_plots || plot_scout_signal) && delta_X == delta_Xs(1)
        f = figure;
        plot(row.signal{1}')
        title('Scouts Signals')
        saveJpegAndFig(f,save_path,char(files_prefix + "- Scouts signal - Delta_T " + delta_T),1)
        close(f)
    end
   
    simulated_recording_sub = simulated_recording(elecrode_nums,:);
        
    if (with_plots || plot_scout_position) && delta_T == delta_Ts(1)
        %plot 2d coordinates - subsample
        f=figure;
        %         text(x_sub, y_sub, string(1:numel(x_sub)));
        text(x_sub, y_sub, z_sub, string(elecrode_nums));
        hold on
        scout_colors = 'rb';
        for scout_num=1:n_scouts
            scout_vertices_inds = row.scouts{1,1}{1,scout_num}.Vertices;
%             scout_vertices_locs = cortex_vars.cortex306716V.Vertices(scout_vertices_inds,:);
            scout_vertices_locs = cortex_vars.cortex.Vertices(scout_vertices_inds,:);
            scatter3(scout_vertices_locs(:,1),scout_vertices_locs(:,2),scout_vertices_locs(:,3),scout_colors(scout_num),'filled');
        end
        axis([min([x_sub; scout_vertices_locs(:,1)]), max([x_sub; scout_vertices_locs(:,1)]), min([y_sub; scout_vertices_locs(:,2)]), max([y_sub; scout_vertices_locs(:,2)]),min([z_sub; scout_vertices_locs(:,3)]), max([z_sub; scout_vertices_locs(:,3)])])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title(["Distance " row.distance])
        saveJpegAndFig(f,save_path,char(files_prefix + "- electode position 3d - subset - Delta_X " + delta_X),1)
        close(f)
    end
    
    En = reshape(1:length(elecrode_nums),size(electrode_map,1),size(electrode_map,2));

%     HT=hilbert(simulated_recording_sub').';
%     HTabs=abs(HT);
%     HTangle=angle(HT);
%     [crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);
   
    max_peaks = 3;
    all_pks = zeros(size(simulated_recording_sub,1),max_peaks);
    all_locs = all_pks;

    for k=1:size(simulated_recording_sub,1)
        [pks,locs] = findpeaks(simulated_recording_sub(k,:));
        all_pks(k,1:length(pks)) = pks;
        all_locs(k,1:length(locs)) = locs;
    end
    
    if percentile
        x_sub = x(elecrode_nums);
        y_sub = y(elecrode_nums);
        z_sub = z(elecrode_nums); 
        all_locs_BU = all_locs;
        min_locs = quantile(all_locs(:,1),percentile);
        max_locs = quantile(all_locs(:,1),1-percentile);
        outliers_inds = find(~(all_locs(:,1)>=max_locs | all_locs(:,1)<=min_locs));
        all_locs = all_locs(outliers_inds,1);
        x_sub = x_sub(outliers_inds);
        y_sub = y_sub(outliers_inds);
        z_sub = z_sub(outliers_inds);
    end
    
    if with_plots || plot_signal
        %plot signal
        f=figure;
        plot(simulated_recording_sub')
        colororder(column_colors(electrode_columns,:))
        lines = findobj(gca,'Type','line');
        first_line_per_group = lines(1:n_electrode_columns:numel(electrode_columns));
        hold on
        scatter(all_locs(:,1),all_pks(:,1),[],column_colors(electrode_columns,:))
        title('Channel subsample signals')
        legend(flipud(first_line_per_group),string(1:4))
        saveJpegAndFig(f,save_path,char(files_prefix + "Signal - first maxima - Delta_T " + delta_T + " Delta_X " + delta_X),1)
        close(f)
    end
    
%     if with_plots
%         %plot signal at half time
%         [pk,loc1] = findpeaks(row.signal{1}(1,:))
%         [pk,loc2] = findpeaks(row.signal{1}(2,:));
%         half_time = round(mean([loc1 loc2]));
%         startEndWave=[0 1];
%         [~,f,~] = plotCrossingsPhysical(simulated_recording_sub(:,half_time),startEndWave,flipud(En),zeros(size(all_pks(:,1)))+10,'Units','frames');
%         end
    
    if with_plots
        %plot PLM
        startEndWave=[1 size(simulated_recording_sub,2)]; %twoGausses
        [~,f,~] = plotCrossingsPhysical(all_locs(:,1),startEndWave,flipud(En),all_pks(:,1),'Units','frames');
        saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map - First Maxima - With Size - Delta_T " + delta_T + " Delta_X " + delta_X),1)
        close(f)
        [~,f,~] = plotCrossingsPhysical(all_locs(:,1),startEndWave,flipud(En),zeros(size(all_pks(:,1)))+30,'Units','frames');
        saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map - First Maxima - No Size - Delta_T " + delta_T + " Delta_X " + delta_X),1)
        
%         [~,f,~] = plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames');
%         saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map With Size- Delta_T " + delta_T + " Delta_X " + delta_X),1)
%         close(f)
%         [~,f,~] = plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),[],'Units','frames');
%         saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map No Size- Delta_T " + delta_T + " Delta_X " + delta_X),1)
        close(f)
    end
    
    if interpolate
        F = scatteredInterpolant([x_sub,y_sub],all_locs(:,1));
        [xq,yq] = meshgrid(linspace(min(x_sub), max(x_sub),20),linspace(min(y_sub), max(y_sub),20));
        vq = F(xq,yq);
        all_locs = vq(:);
        f=figure;
        surf(xq,yq,vq)
        saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map - Interpolated surf - Delta_T " + delta_T + " Delta_X " + delta_X),1)
        close(f)
    end
    
    phase_space_row = find(delta_Ts==delta_T);
    phase_space_col = find(delta_Xs==delta_X);
    
    %PLDC
    if interpolate
        cors(phase_space_row,phase_space_col) = corr(xq(:),all_locs(:,1));
    else
        cors(phase_space_row,phase_space_col) = corr(x_sub,all_locs(:,1));
    end
    
    [dipLFP, dip_p_values(phase_space_row,phase_space_col)] = hartigansdipsigniftest(sort(all_locs(:,1)), 500);
%         cors(i,j) = corr(x_sub,crossings{1}(:,1));
%         [dipLFP, dip_p_values(i,j)] = hartigansdipsigniftest(sort(crossings{1}(:,1)), 500);
    t_sanity(phase_space_row,phase_space_col) = delta_T;
    x_sanity(phase_space_row,phase_space_col) = delta_X;
    
    if with_plots
        f=figure;
%         scatter(x_sub,crossings{1}(:,1),'filled')
        if interpolate
            scatter(xq(:),all_locs(:,1),'filled')
        else
            scatter(x_sub,all_locs(:,1),'filled')
        end
        xlabel('x')
        ylabel('t')
        title(['PLDC - ' num2str(cors(phase_space_row,phase_space_col))])
        saveJpegAndFig(f,save_path,char("PLDC - first maxima - Delta_T " + delta_T + " Delta_X " + delta_X),1)
        close(f)

        %dip
        f=figure;

        histogram(all_locs(:,1),20)
        title(['DIP p-value: ' num2str(dip_p_values(phase_space_row,phase_space_col))])
        saveJpegAndFig(f,save_path,char("crossing times hist and dip - first maxima - Delta_T " + delta_T + " Delta_X " + delta_X'),1)
        close(f)
    end   
end

