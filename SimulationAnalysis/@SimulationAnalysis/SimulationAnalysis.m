classdef SimulationAnalysis
    %SIMULATIONRUN Class for handling a single simulation analysis.
    %   Contains all relevant simulation data (cortex, sources, etc.),
    %   spatiotemporal setting (delta_X,delta_T) and result metrics
    %   (PLDC,dip_p_value)
    
    properties
        cortex %cortex struct including cortex name and vertices_positions
        signal_source %signal source object - scouts or dipoles
        electrode_set % electrode set object - containing positions and recorded signals
        electrode_grid % aka En - 2d grid of channel numbers
        electrode_subset %usually electrode_grid(:) unless outliers are removed (see percentile plotting param)
        sampling_rate %Simulation sampling_rate in Hz - applies both for signal source and elecrode recorded signal
        source_delta_T %delta_T between sources charactarizing this analysis
        source_delta_X %delta_X between sources charactarizing this analysis
        dip_p_value %the dip-test p-value of this simulation
        PLDC %the phase latency-distance correlation of this simulation
    end
    
    methods
        function obj = SimulationAnalysis(signal_source, electrode_set, electrode_grid, sampling_rate, cortex, source_delta_T, source_delta_X)
            %SIMULATIONRUN Construct an instance of this class
            %   INPUT:
            %       - signal_source -
            %            SignalSource object (Scouts or Dipoles)
            %            used in simulation
            %       - electrode_signals [n_electrodes X n_samples]
            %            the signals recorded by the electrodes
            %       - electrode_positions - n_electrodex X 3
            %       - electrode_grid - aka En
            %       - source_delta_T, source_delta_X - characteristics
            %            identifying this simulation's signal source
            
            obj.signal_source = signal_source;
            obj.electrode_set = electrode_set; %ElectrodeSet(electrode_signals, electrode_positions, electrode_grid);
            obj.electrode_grid = electrode_grid;
            obj.electrode_subset = obj.electrode_grid(:);
            obj.sampling_rate = sampling_rate;
            obj.cortex = cortex;
            obj.source_delta_T = source_delta_T;
            obj.source_delta_X = source_delta_X;
        end

        function [phase_event_heights, phase_event_locs] = get_phase_events(obj,event_type,start_end_wave)
            % Currently only event_type="first_max" is supported
            if event_type~="first_max"
                error('Currently only event_type="first_max" is supported')
            end

            max_peaks = 3;
            all_pks = zeros(length(obj.electrode_subset),max_peaks);
            all_locs = all_pks;

            for k=1:length(obj.electrode_subset)
                signal_in_wave = obj.electrode_set.electrode_signals(obj.electrode_subset(k),start_end_wave(1):start_end_wave(2)); 
                [pks,locs] = findpeaks(signal_in_wave);
                all_pks(k,1:length(pks)) = pks;
                all_locs(k,1:length(locs)) = locs + start_end_wave(1) - 1;
            end

            phase_event_heights = all_pks(:,1);
            phase_event_locs = all_locs(:,1);
        end
        
        function [projected] = project_2d_to_1d(obj,axis_1,axis_2)
            [~,score,~] = pca([axis_1, axis_2]);
            projected = score(:,1);
        end
        function [phase_event_heights, phase_event_locs, new_electrode_subset] = remove_events_outliers(obj,phase_event_heights,phase_event_locs,percentile)
                min_locs = quantile(phase_event_locs,percentile);
                max_locs = quantile(phase_event_locs,1-percentile);
                outliers_inds = find(~(phase_event_locs>=max_locs | phase_event_locs<=min_locs));
                phase_event_locs = phase_event_locs(outliers_inds);
                phase_event_heights = phase_event_heights(outliers_inds);
                new_electrode_subset = obj.electrode_subset(outliers_inds);
        end

        function dip_p_value = calc_dip_p_value(obj, phase_event_locs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [~, dip_p_value] = hartigansdipsigniftest(sort(phase_event_locs), 500);
        end
        
        function PLDC = calc_PLDC(obj, phase_event_locs)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            PLDC = corr(obj.electrode_set.electrode_positions(obj.electrode_subset,1),phase_event_locs);
        end

        [dip_p_value, PLDC] = run_single_analysis(obj, plotting_config, save_path, files_prefix, is_first_temporal, is_first_spatial)

        function En = electrode_grid_to_En(obj)
            subset_mask = ismember(obj.electrode_grid,obj.electrode_subset);
            En = NaN(size(obj.electrode_grid));
            En(subset_mask) = 1:sum(subset_mask,'all');
        end
    end
end

