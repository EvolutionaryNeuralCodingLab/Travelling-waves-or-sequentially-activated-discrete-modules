function [dip_p_value, PLDC] = run_single_analysis(obj, plotting_config, save_path, files_prefix, is_first_temporal, is_first_spatial)
%ANALYZE_SINGLE_SIMULATION Summary of this function goes here
%   Detailed explanation goes here

% Plotting config with_plots, plot_summary, plot_signal, plot_source_signal, plot_source_position, percentile, interpolate

% Get Phase events phase_event_locs, 

%Fill missing plotting values with defaults
plotting_config = fill_missing_configs_with_defaults(plotting_config,get_default_plotting_config());

if isempty(plotting_config.start_end_wave)
    analysis_start_end_wave = [1 size(obj.electrode_set.electrode_signals,2)];
else
    analysis_start_end_wave = plotting_config.start_end_wave;
end

[phase_event_heights, phase_event_locs] = obj.get_phase_events("first_max",analysis_start_end_wave);
if plotting_config.percentile
    [phase_event_heights, phase_event_locs, obj.electrode_subset] = obj.remove_events_outliers(phase_event_heights,phase_event_locs,plotting_config.percentile);
end

En = obj.electrode_grid_to_En();

%plot 3d coordinates all electrodes
if (plotting_config.with_plots || plotting_config.plot_source_position) && is_first_temporal
    f=figure;
    ax = obj.electrode_set.plot_electrode_positions(gca);
    hold on
    ax = obj.signal_source.plot_source_positions(ax);
    [electrode_x_lims, electrode_y_lims, electrode_z_lims] = obj.electrode_set.get_min_max_positions();
    axis([min(obj.signal_source.x_lims(1), electrode_x_lims(1)), ...
          max(obj.signal_source.x_lims(2), electrode_x_lims(2)), ...
          min(obj.signal_source.y_lims(1), electrode_y_lims(1)), ...
          max(obj.signal_source.y_lims(2), electrode_y_lims(2)), ...
          min(obj.signal_source.z_lims(1), electrode_z_lims(1)), ...
          max(obj.signal_source.z_lims(2), electrode_z_lims(2))])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['Elecrode Position (', obj.signal_source.SOURCE_TYPE, ' in color) - Delta X ', num2str(obj.source_delta_X)])
    saveJpegAndFig(f,save_path,char(files_prefix + "electode position 3d - Delta X " + obj.source_delta_X),1)
    close(f)
end

if (plotting_config.with_plots || plotting_config.plot_source_signal) && is_first_spatial
    f = figure;
    ax = obj.signal_source.plot_source_signal(plotting_config.start_end_wave, gca); %plotting_config.start_end_wave sent and not analysis_start_end_wave so it'll be empty if needed
%     obj.convert_x_from_sample_to_seconds(ax);
    title(obj.signal_source.SOURCE_TYPE + " Signals - Delta T " + obj.source_delta_T)
    xlabel("Time [s]")
    saveJpegAndFig(f,save_path,char(files_prefix + obj.signal_source.SOURCE_TYPE + " signal - Delta_T " + obj.source_delta_T),1)
    close(f)
end

if (plotting_config.with_plots || plotting_config.plot_source_position) && is_first_temporal
    f=figure;
    ax = obj.electrode_set.plot_electrode_positions(gca,obj.electrode_grid(:));
    hold on
    ax = obj.signal_source.plot_source_positions(ax);
    [electrode_x_lims, electrode_y_lims, electrode_z_lims] = obj.electrode_set.get_min_max_positions(obj.electrode_grid(:));
    axis([min(obj.signal_source.x_lims(1), electrode_x_lims(1)), ...
          max(obj.signal_source.x_lims(2), electrode_x_lims(2)), ...
          min(obj.signal_source.y_lims(1), electrode_y_lims(1)), ...
          max(obj.signal_source.y_lims(2), electrode_y_lims(2)), ...
          min(obj.signal_source.z_lims(1), electrode_z_lims(1)), ...
          max(obj.signal_source.z_lims(2), electrode_z_lims(2))])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['Elecrode Position (', obj.signal_source.SOURCE_TYPE, ' in color) - Delta X ', num2str(obj.source_delta_X)])
    saveJpegAndFig(f,save_path,char(files_prefix + "electode position 3d - subset - Delta X " + num2str(obj.source_delta_X)),1)
    close(f)
end

if plotting_config.with_plots || plotting_config.plot_signal
    %plot signal
    f=figure;
    ax = obj.electrode_set.plot_signal(gca, obj.electrode_grid, obj.electrode_subset, plotting_config.start_end_wave, phase_event_locs, phase_event_heights); %plotting_config.start_end_wave sent and not analysis_start_end_wave so it'll be empty if needed
    xlabel('Time [s]')
    title("Signal - first maxima - Delta T " + obj.source_delta_T + " Delta X " + obj.source_delta_X)
    saveJpegAndFig(f,save_path,char(files_prefix + "Signal - first maxima - Delta_T " + obj.source_delta_T + " Delta_X " + obj.source_delta_X),1)
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

%%TODO: implement via objects + fix colorbar
if plotting_config.with_plots
    %plot PLM
%     startEndWave=[1 size(obj.electrode_set.electrode_signals,2)]; %twoGausses
    [~,f,~] = plotCrossingsPhysical(phase_event_locs * 1000 / obj.sampling_rate,analysis_start_end_wave * 1000 / obj.sampling_rate,flipud(En),phase_event_heights,'Units','ms');
    title("Phase Latency Map - First Maxima - Delta T " + obj.source_delta_T + " Delta X " + obj.source_delta_X)
    saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map - First Maxima - With Size - Delta_T " + obj.source_delta_T + " Delta_X " + obj.source_delta_X),1)
    close(f)
    [~,f,~] = plotCrossingsPhysical(phase_event_locs * 1000 / obj.sampling_rate,analysis_start_end_wave * 1000 / obj.sampling_rate,flipud(En),zeros(size(phase_event_heights))+30,'Units','ms');
    title("Phase Latency Map - First Maxima - Delta T " + obj.source_delta_T + " Delta X " + obj.source_delta_X)
    saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map - First Maxima - No Size - Delta_T " + obj.source_delta_T + " Delta_X " + obj.source_delta_X),1)
    
    %         [~,f,~] = plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),hilbertAmps{1},'Units','frames');
    %         saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map With Size- Delta_T " + delta_T + " Delta_X " + delta_X),1)
    %         close(f)
    %         [~,f,~] = plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),[],'Units','frames');
    %         saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map No Size- Delta_T " + delta_T + " Delta_X " + delta_X),1)
    close(f)
end

%%TODO
% if plotting_config.interpolate
%     F = scatteredInterpolant([x_sub,y_sub],all_locs(:,1));
%     [xq,yq] = meshgrid(linspace(min(x_sub), max(x_sub),20),linspace(min(y_sub), max(y_sub),20));
%     vq = F(xq,yq);
%     all_locs = vq(:);
%     f=figure;
%     surf(xq,yq,vq)
%     saveJpegAndFig(f,save_path,char(files_prefix + "Phase Latency Map - Interpolated surf - Delta_T " + delta_T + " Delta_X " + delta_X),1)
%     close(f)
% end


dip_p_value = obj.calc_dip_p_value(phase_event_locs);
PLDC = obj.calc_PLDC(phase_event_locs);

% %PLDC
% if interpolate
%     cors(phase_space_row,phase_space_col) = corr(xq(:),phase_event_locs(:,1));
% else
%     
% end

% [dipLFP, dip_p_values(phase_space_row,phase_space_col)] = hartigansdipsigniftest(sort(phase_event_locs(:,1)), 500);
%         cors(i,j) = corr(x_sub,crossings{1}(:,1));
%         [dipLFP, dip_p_values(i,j)] = hartigansdipsigniftest(sort(crossings{1}(:,1)), 500);


if plotting_config.with_plots
    f=figure;
    %         scatter(x_sub,crossings{1}(:,1),'filled')
%     if interpolate
%         scatter(xq(:),phase_event_locs(:,1),'filled')
%     else
    scatter(obj.electrode_set.electrode_positions(obj.electrode_subset,1),phase_event_locs / obj.sampling_rate,'filled')
%     end
    xlabel('X')
    ylabel('Time [s]')
    title(char("PLDC " + num2str(PLDC) + " - first maxima - Delta T " + obj.source_delta_T + " Delta X " + obj.source_delta_X))
    saveJpegAndFig(f,save_path,char("PLDC - first maxima - Delta_T " + obj.source_delta_T + " Delta_X " + obj.source_delta_X),1)
    close(f)
    
    %dip
    f=figure;
    
    histogram(phase_event_locs / obj.sampling_rate,20)
    xlabel('Phase Latencies [s]')
    ylabel('Frequency')
    title(char("DIP p-value: " + num2str(dip_p_value) + " Delta T " + obj.source_delta_T + " Delta X " + obj.source_delta_X'))
    saveJpegAndFig(f,save_path,char("crossing times hist and dip - first maxima - Delta_T " + obj.source_delta_T + " Delta_X " + obj.source_delta_X'),1)
    close(f)
end
end

