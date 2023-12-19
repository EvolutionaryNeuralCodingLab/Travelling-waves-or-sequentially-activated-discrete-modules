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

if plotting_config.duplicate_hemisphere
   EEGcap = load('/media/sil2/Data/Yuval O/Rotem Simulations/SimulationAnalysis/spatial_data/EEGcap.mat'); %TODO: Move path to consts 
   EEGcap = EEGcap.EEGcap;
   [simulation_table, ~] = find_electrode_couples(EEGcap, simulation_table);
end

for i=1:n_rows
   row = simulation_table(i,:);
   delta_T = row.deltaT;
   delta_X = row.distance;
   electrode_signals = row.EEG_recordings{1};
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

   is_first_temporal = delta_T == delta_Ts(1);
   is_first_spatial = delta_X == delta_Xs(1);

   phase_space_row = find(delta_Ts==delta_T);
   phase_space_col = find(delta_Xs==delta_X);
   [dip_p_values(phase_space_row,phase_space_col), PLDCs(phase_space_row,phase_space_col)] = single_run.run_single_analysis(plotting_config, save_path, files_prefix, is_first_temporal, is_first_spatial);
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
end