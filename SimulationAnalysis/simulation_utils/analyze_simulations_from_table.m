function [dip_p_values, PLDCs, phase_space_dT, phase_space_dX, delta_Ts, delta_Xs] = analyze_simulations_from_table(simulation_table, electrode_channel_poitions, electrode_grid, source_type, sampling_rate, cortex, plotting_config, save_path, files_prefix)
%ANALYZE_SIMULATIONS Analyzes all simulations save in a table format.
%   Both returns and saves metrics! Unless save_path='' (is_empty) then it
%   doesn't save.


%Fill missing plotting values with defaults
plotting_config = fill_missing_configs_with_defaults(plotting_config,get_default_plotting_config());

delta_Ts = unique(simulation_table.deltaT);
delta_Xs = unique(simulation_table.distance);

n_rows = height(simulation_table);

PLDCs = zeros(length(delta_Ts),length(delta_Xs));
dip_p_values = PLDCs;
phase_space_dT = PLDCs;
phase_space_dX = PLDCs;

% elecrode_nums = electrode_map(~isnan(electrode_map));
% nCh_sub = length(elecrode_nums);
% 
% electrode_columns = repmat(1:size(electrode_map,2),size(electrode_map,1),1);
% n_electrode_columns = size(electrode_map,2);
% column_colors = parula(n_electrode_columns);

% % channel coordinates
% x = channel_locations(:,1);
% y = channel_locations(:,2);
% z = channel_locations(:,3);        
% x_sub = x(elecrode_nums);
% y_sub = y(elecrode_nums);
% z_sub = z(elecrode_nums);

if plotting_config.duplicate_hemisphere
   EEGcap = load('/media/sil2/Data/Yuval O/Rotem Simulations/SimulationAnalysis/spatial_data/EEGcap.mat'); %TODO: Move path to consts 
   EEGcap = EEGcap.EEGcap;
   [simulation_table, ~] = find_electrode_couples(EEGcap, simulation_table); %TODO: Refactor Rotem's code
end

for i=1:n_rows
   row = simulation_table(i,:);
   delta_T = row.deltaT;
   delta_X = row.distance;
%    disp(['delta_T ' num2str(delta_T) ' delta_X ' num2str(delta_X)])
   electrode_signals = row.EEG_recordings{1};
    %Above is Rotems' code to flip all channels
    %This is old code for flipping only En electrode. TODO: Delete after starting code versioning   
%    if plotting_config.duplicate_hemisphere
%        temp_signal = row.EEG_recordings{1};
%        flipped_electrode_grid = flipud(electrode_grid);
%        temp_signal_flipped = temp_signal;
%        temp_signal_flipped(electrode_grid(:),:) = temp_signal_flipped(flipped_electrode_grid(:),:);
%        temp_signal_duplicated = temp_signal;
%        temp_signal_duplicated(electrode_grid(:),:) = temp_signal(electrode_grid(:),:) + temp_signal_flipped(electrode_grid(:),:);
%        electrode_signals = temp_signal_duplicated;
%    else
%        electrode_signals = row.EEG_recordings{1};
%    end
   crop_at = get_crop_at_if_empty(plotting_config.crop_at,size(electrode_signals,2));
   electrode_set = ElectrodeSet(electrode_signals(:,crop_at(1):crop_at(2)), electrode_channel_poitions,plotting_config);

   if plotting_config.ignore_source
       sources = [];
   else
       if strcmp(source_type,"scouts")
        sources = Scouts.init_source_from_table_row(row, cortex);
       elseif strcmp(source_type,"dipoles")
        sources = Dipoles.init_source_from_table_row(row, cortex);
       else
           error("Unkown source type.")
       end
   end
    
   single_run = SimulationAnalysis(sources, electrode_set, electrode_grid, sampling_rate, cortex, delta_T, delta_X);

   single_run.electrode_set.parent_analysis = single_run;
   single_run.signal_source.parent_analysis = single_run;

%    source_position_col = 'scouts';
%    source_signal = simulation_table(1,"signal").signal{1};
%    vertices_inds = cellfun(@(x) x.Vertices, simulation_table(1,source_position_col).scouts{1}, 'UniformOutput', false);
%    scout_positions = Scouts.
%    Scouts.scouts_positions_from_inds(vertices_inds,cortex.vertices_positions);

%    sources = Scouts(scout_positions,source_signal);

   is_first_temporal = delta_T == delta_Ts(1);
   is_first_spatial = delta_X == delta_Xs(1);

   phase_space_row = find(delta_Ts==delta_T);
   phase_space_col = find(delta_Xs==delta_X);
   [dip_p_values(phase_space_row,phase_space_col), PLDCs(phase_space_row,phase_space_col)] = single_run.run_single_analysis(plotting_config, save_path, files_prefix, is_first_temporal, is_first_spatial);
   %    analyze_single_simulation(signal_source, electrode_set, electrode_grid_for_analysis, source_delta_T, source_delta_X, plotting_config, is_first_temporal, is_first_spatial)
   if plotting_config.save_analysis_objects
       object_save_path = [save_path 'simulation_analysis_objects/'];
       if ~isfolder(object_save_path)
           mkdir(object_save_path)
       end
       save([object_save_path 'simulation_analysis_delta_T_' num2str(delta_T) '_delta_X_' num2str(delta_X) '.mat'],'single_run')
   end
   phase_space_dT(phase_space_row,phase_space_col) = delta_T;
   phase_space_dX(phase_space_row,phase_space_col) = delta_X;
end

%%TODO: Move this to external plotting funcitons
if ~isempty(save_path)
    save(char(save_path + files_prefix + 'wave metrics - first maxima'),'PLDCs','dip_p_values','delta_Ts','delta_Xs','electrode_grid', 'sources')
end


if plotting_config.with_plots || plotting_config.plot_summary
    T_start_end = [1,length(delta_Ts)];
    X_start_end = [1,length(delta_Xs)];
    f = plot_phase_space(dip_p_values,delta_Ts, delta_Xs, T_start_end,X_start_end,1);
    
    if plotting_config.with_plots
        saveJpegAndFig(f,save_path,'phase space DIP dipols - first maxima',1)
        close(f)
    end
    T_start_end = [1,length(delta_Ts)];
    X_start_end = [1,length(delta_Xs)];
    if isfield(plotting_config,'PLDC_summary_scale')
        img_range = plotting_config.PLDC_summary_scale;
    else
        img_range = [0.8,1];
    end
    f = plot_phase_space(PLDCs,delta_Ts, delta_Xs, T_start_end,X_start_end,0,img_range);
    if plotting_config.with_plots
        saveJpegAndFig(f,save_path,'phase space PLDC dipols - first maxima - colorbar adjusted',1)
        close(f)
    end
        
        
%         T_start_end = [2,length(delta_Ts)];
%         X_start_end = [3,length(delta_Xs)];
% 
%         f=figure;    
%         imagesc(dip_p_values(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2)))
%         colorbar()
%         yticks(1:(T_start_end(2)-T_start_end(1)+1))
%         yticklabels(num2cell(delta_Ts(T_start_end(1):T_start_end(2))))
%         ylabel('Delta T / std(T)')
%         xticks(1:(X_start_end(2)-X_start_end(1)+1))
%         xticklabels(num2cell(delta_Xs(X_start_end(1):X_start_end(2))))
%         xlabel('dipole dist / electrode dist')
%         set(gca,'YDir','normal')
%         hold on
%         [row,col] = ind2sub(size(dip_p_values(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))),find(dip_p_values(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))<0.1));
%         scatter(col,row,'k','filled')
%         [row,col] = ind2sub(size(dip_p_values(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))),find(dip_p_values(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))<0.05));
%         scatter(col,row,'white')
%         saveJpegAndFig(f,save_path,'phase space DIP dipols - first maxima - partial',1)
 end