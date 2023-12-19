classdef SignalSource
    %SIGNALSOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Abstract, Constant)
        SOURCE_TYPE
    end
    
    properties (Abstract)
        n_sources
        sources_positions %cell array with n_signal_sources cells for each source
        x_lims
        y_lims
        z_lims
        signal
    end

    properties
        parent_analysis %Will be assigned once exists. sampling_rate comes from here
    end
    
    methods (Abstract)
        ax = plot_source_positions(obj, ax)
        %Given an ax handle ax, plots the location of source. Returns
        %handle
        [x, y, z] = get_source_positions(obj, source_num)
        %Return the x,y,z position (or positions, if scout) for the
        %source_num's source object (source_num<=n_signal_sources).
        % length(x)=length(y)=length(z)=1 for dipoles, or the number of
        % vertices in a scout for scouts.
        
    end
    
    methods (Abstract, Static)
        [source] = init_source_from_table_row(obj, row);
    end

    methods
        function ax = plot_source_signal(obj, start_end_wave, ax)
            t = (1:size(obj.signal,2)) / obj.parent_analysis.sampling_rate;
            plot(ax, t, obj.signal')
            if ~isempty(start_end_wave)
                hold on
                xline(start_end_wave(1) / obj.parent_analysis.sampling_rate,'--k');
                xline(start_end_wave(2) / obj.parent_analysis.sampling_rate,'--k');
                legend({'Source signal 1', 'Source signal 2', 'Wave Start', 'Wave End'})
            else
                legend({'Source signal 1', 'Source signal 2'})
            end
        end
    end
end

