classdef ElectrodeSet
    %ELECTRODESET Object for handling simulated recording electrodes.
    
    properties
        n_electrodes
        electrode_positions
        electrode_signals
        electrode_nums
        raw_recorded_signal %exists only if ws supplied
        recorded_signal_with_noise %exists only if noise_gauss_params supplied
        filter_obj %exists only if ws supplied
        parent_analysis %Will be assigned once exists. sampling_rate comes from here
    end
    
    methods
        function obj = ElectrodeSet(electrode_signals, electrode_positions, plotting_config)
            %ELECTRODESET Construct an instance of this class
            %   electrode_grid aka En
            % If ws (1x2) is passed, electrode_signals are bandpass filtered, keeping
            % the raw recorded signal in obj.raw_recorded_signal and the
            % filter object in obj.filter_obj
            % If fs (sampling frequency, Hz) also passed, ws will be
            % treated as Hz
            % If pad_samples is provided, it duplicates the start and end
            % pad_samples samples n_dupes times
            
            %Fill missing plotting values with defaults
            plotting_config = fill_missing_configs_with_defaults(plotting_config,get_default_plotting_config());
            if ~isempty(plotting_config.noise_gauss_params)
                obj.recorded_signal_with_noise = electrode_signals .* ...
                                                 normrnd(plotting_config.noise_gauss_params(1),...
                                                         plotting_config.noise_gauss_params(2),...
                                                         size(electrode_signals));
            end
            if ~isempty(plotting_config.ws)
                if ~isempty(plotting_config.noise_gauss_params)
                    signal_to_filter = obj.recorded_signal_with_noise;
                else
                    signal_to_filter = electrode_signals;
                end
                %Prepare for filtering
                padding_block_start = signal_to_filter(:,1:plotting_config.pad_samples);
                padding_block_end = signal_to_filter(:,(size(signal_to_filter,2)-plotting_config.pad_samples+1):end);
                if plotting_config.flip_signal
                    padding_block_start = fliplr(padding_block_start);
                    padding_block_end = fliplr(padding_block_end);
                end
                temp_signals = [repmat(padding_block_start,1,plotting_config.n_reps) signal_to_filter repmat(padding_block_end,1,plotting_config.n_reps)];
                %Filtering the data. First check if cache exists
                cfs = CachedFilteredSignal();
%                 filtered_signal_cache_struct = cfs.create_filtered_signal_cache_struct(electrode_signals,plotting_config.ws, plotting_config.fs, plotting_config.crop_at, plotting_config.pad_samples, plotting_config.flip_signal, plotting_config.n_reps);
                filtered_signal_cache_struct = cfs.create_filtered_signal_cache_struct(temp_signals,plotting_config.ws, plotting_config.fs);
                hached_full_path = cfs.get_hached_full_path(filtered_signal_cache_struct);
                if isfile(hached_full_path)
                    cached_filtered = load(hached_full_path);
                    filtered_data = cached_filtered.filtered_data;
                    filter_obj = cached_filtered.filter_obj;
                else
                    disp('Filtering signal...')
                    if ~isempty(plotting_config.fs)
                        [filtered_data, filter_obj] = bandpass(temp_signals',plotting_config.ws,plotting_config.fs);
                    else
                        [filtered_data, filter_obj] = bandpass(temp_signals',plotting_config.ws);
                    end
                    disp('Done filtering...')
                    disp(['Saving in ' hached_full_path])
                    if ~exist(cfs.cached_folder, 'dir')
                       mkdir(cfs.cached_folder)
                    end
                    save(hached_full_path,'filtered_data', 'filter_obj')
                end
                filtered_data = filtered_data((plotting_config.pad_samples*plotting_config.n_reps+1):(end-(plotting_config.pad_samples*plotting_config.n_reps)),:);
                obj.electrode_signals = filtered_data';
                obj.raw_recorded_signal = electrode_signals;
                obj.filter_obj = filter_obj;
            else
                if ~isempty(plotting_config.noise_gauss_params)
                    obj.raw_recorded_signal = electrode_signals;
                    obj.electrode_signals =  obj.recorded_signal_with_noise;
                else
                    obj.electrode_signals = electrode_signals;
                end
            end
            obj.electrode_positions = electrode_positions;
            obj.n_electrodes = size(electrode_positions,1);
            obj.electrode_nums = 1:obj.n_electrodes;
%             if exist(electrode_grid,'var')
%                 obj.electrode_grid = electrode_grid;
%                 obj.electrode_nums = electrode_grid(:);
%                 obj.has_grid = True;
%             else
%                 obj.electrode_grid = [];
%                 obj.electrode_nums = 1:n_electrodes;
%                 obj.has_grid = False;
%             end
        end
        
        function ax = plot_electrode_positions(obj, ax, subset)
            if ~exist('subset','var')
                subset = obj.electrode_nums;
            end
            x = obj.electrode_positions(subset,1);
            y = obj.electrode_positions(subset,2);
            z = obj.electrode_positions(subset,3);

            text(ax, x,y,z, string(subset(:)));
            ax = gca;
        end
        
        function [ax] = plot_signal(obj, ax, electrode_grid, subset, start_end_wave, phase_event_locs, phase_event_heights)
            %METHOD1 Summary of this method goes here
            %   TODO: allow subset to be empty/not delievered

            if ~exist('ax','var')
                figure;
                ax=gca;
            end
            if isempty(obj.parent_analysis)
                sampling_rate = 1; %#Todo: get this as optional varain
            else
                sampling_rate = obj.parent_analysis.sampling_rate;
            end
            t = (1:size(obj.electrode_signals,2)) / sampling_rate;
            if ~exist('electrode_grid','var')
                ax = plot(t,obj.electrode_signals');
            else
                electrode_columns = repmat(1:size(electrode_grid,2),size(electrode_grid,1),1);
                n_electrode_columns = size(electrode_grid,2);
                column_colors = parula(n_electrode_columns);
                subset_mask = ismember(electrode_grid,subset);

                plot(t, obj.electrode_signals(subset,:)')
                colororder(column_colors(electrode_columns(subset_mask),:))
                lines = findobj(gca,'Type','line');
                first_lines_nums = cumsum([1 sum(subset_mask,1)]);
                first_line_per_group = lines(first_lines_nums(1:end-1));
                if exist('phase_event_locs','var')  && exist('phase_event_heights','var')
                    hold on
                    scatter(phase_event_locs / obj.parent_analysis.sampling_rate,phase_event_heights,[],column_colors(electrode_columns(subset_mask),:))
                end
                legend(flipud(first_line_per_group),string(1:n_electrode_columns))
            end
            if exist('start_end_wave','var') && ~isempty(start_end_wave)
                hold on
                xline(start_end_wave(1) / obj.parent_analysis.sampling_rate,'--k','DisplayName','Wave Start');
                xline(start_end_wave(2) / obj.parent_analysis.sampling_rate,'--k','DisplayName','Wave End');
            end
        end

        function [x_lims, y_lims, z_lims] = get_min_max_positions(obj, subset)
            if ~exist('subset','var')
                subset = obj.electrode_nums;
            end
            x = obj.electrode_positions(subset,1);
            y = obj.electrode_positions(subset,2);
            z = obj.electrode_positions(subset,3);
            x_lims = [min(x), max(x)];
            y_lims = [min(y), max(y)];
            z_lims = [min(z), max(z)];
        end
    end
end

