classdef Scouts < SignalSource
    %SCOUTS Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        SOURCE_TYPE = 'Scouts'
        SCOUT_COLORS = 'brgyk'
    end
    properties
        n_sources
        sources_positions %cell array with n_signal_sources cells for each source
        x_lims
        y_lims
        z_lims
        signal
    end
       
    methods
        function obj = Scouts(scouts_positions, signal)
            %SIGNALSOURCE Construct an instance of this class
            %   scouts_positions: cell array containing the x,y,z positions of the scouts.
            %   signal: array containing the scouts signal (n_sourcesXn_samples)
            obj.n_sources = length(scouts_positions);
            obj.sources_positions = scouts_positions;
            obj.signal = signal;
            
            obj.x_lims = [inf -inf];
            obj.y_lims = [inf -inf];
            obj.z_lims = [inf -inf];
            for i=1:obj.n_sources
                [x,y,z] = obj.get_source_positions(i);
                obj.x_lims = [min(obj.x_lims(1), min(x)), max(obj.x_lims(2), max(x))];
                obj.y_lims = [min(obj.y_lims(1), min(y)), max(obj.y_lims(2), max(y))];
                obj.z_lims = [min(obj.z_lims(1), min(z)), max(obj.z_lims(2), max(z))];
            end
        end

        function [x, y, z] = get_source_positions(obj, source_num)
            position_array = obj.sources_positions{source_num};
            [x, y, z] = deal(position_array(:,1), position_array(:,2), position_array(:,3));
        end
        
        function ax = plot_source_positions(obj, ax, flip_y, flip_around)
            if ~exist('flip_y','var')
                flip_y = false;
            end
            if ~exist('flip_around','var')
                flip_around = 0;
            end
            hold(ax, 'on')
            if flip_y
                for scout_num=1:obj.n_sources
%                     %gray
%                     scatter3(ax, obj.sources_positions{scout_num}(:,1),-obj.sources_positions{scout_num}(:,2),obj.sources_positions{scout_num}(:,3),[],repmat([0.5 0.5 0.5],numel(obj.sources_positions{scout_num}(:,1)),1));
%                     %gray filled
%                     scatter3(ax, obj.sources_positions{scout_num}(:,1),-obj.sources_positions{scout_num}(:,2),obj.sources_positions{scout_num}(:,3),[],repmat([0.5 0.5 0.5],numel(obj.sources_positions{scout_num}(:,1)),1),'filled');
                    %colored filled
                    scatter3(ax, obj.sources_positions{scout_num}(:,1),2*flip_around - obj.sources_positions{scout_num}(:,2),obj.sources_positions{scout_num}(:,3),obj.SCOUT_COLORS(scout_num),'filled');
%                     %colored hollow
%                     scatter3(ax, obj.sources_positions{scout_num}(:,1),-obj.sources_positions{scout_num}(:,2),obj.sources_positions{scout_num}(:,3),obj.SCOUT_COLORS(scout_num));
                end
            else
                for scout_num=1:obj.n_sources
                    scatter3(ax, obj.sources_positions{scout_num}(:,1),obj.sources_positions{scout_num}(:,2),obj.sources_positions{scout_num}(:,3),obj.SCOUT_COLORS(scout_num),'filled');
                end
            end
        end
    end

    methods(Static)
        function [source] = init_source_from_table_row(row, cortex)
            n_scouts = length(row.scouts{1});
            scouts_positions = cell(n_scouts,1);
            for scout_num=1:n_scouts
                scouts_vertices_inds = row.scouts{1,1}{1,scout_num}.Vertices;
                scouts_positions{scout_num} = cortex.vertices_positions(scouts_vertices_inds,:);
%                 scatter3(scout_vertices_locs(:,1),scout_vertices_locs(:,2),scout_vertices_locs(:,3),scout_colors(scout_num),'filled');
                signal = row.signal{1};
            end
            source = Scouts(scouts_positions, signal);
        end
       
       function scouts_positions = scouts_positions_from_inds(vertices_inds, vertices_positions)
       %   vertices_inds: cell array of scouts indices (cell per scout).
       %   vertices_positions: vertices_positions(vertices_inds{1}(1),:) are
       %   the x,y,z positions of the first scout
           n_scouts = length(vertices_inds);
           scouts_positions = cell(1,n_scouts);
           for i=1:n_scouts
               scout_inds = vertices_inds{i};
               scouts_positions{i} = vertices_positions(scout_inds,:);
           end
       end
    end    

end

