function [all_scouts_dipoles_positions, mean_electrode_dist, scouts_dist_per_electrode_dist, scouts_stds, scouts_dist_per_std] = ...
    calc_distances(dxs, dt, electrode_grid, cortex, simulation_table, electrode_channel_poitions)
%CALC_DISTANCES calculates the distance of scouts.
%   dt varargin is arbitrary, just needs to be present in table.

all_scouts_dipoles_positions = cell(length(dxs),1);
scouts_mean_distances = zeros(length(dxs),1);
scouts_stds = zeros(length(dxs),1);

% Get scouts positions and distances
for i = 1:length(dxs)
    dx = dxs(i);
    %%% Get relevant row
    row_num = get_row_num_by_dt_dx(simulation_table,dt,dx);
    row = simulation_table(row_num ,:);
    sources = Scouts.init_source_from_table_row(row, cortex);
    
    %save positions
    all_scouts_dipoles_positions{i} = sources.sources_positions;
    
    % Mean scouts position
    scout1_mean_pos = mean(sources.sources_positions{1, 1});
    scout1_std = std(sources.sources_positions{1, 1}(:,1));
    scout2_mean_pos = mean(sources.sources_positions{2, 1});
    scout2_std = std(sources.sources_positions{2, 1}(:,1));

    scouts_mean_distances(i) = sqrt(sum((scout2_mean_pos - scout1_mean_pos).^2));
    scouts_stds(i) = mean([scout1_std,scout2_std]);
end

% Get electrode mean distance
electrode_grid_shifted = electrode_grid(:,2:end);
electrode_grid_unshifted = electrode_grid(:,1:end-1);

x1 = electrode_channel_poitions(electrode_grid_unshifted(:), 1);
x2 = electrode_channel_poitions(electrode_grid_shifted(:), 1);
y1 = electrode_channel_poitions(electrode_grid_unshifted(:), 2);
y2 = electrode_channel_poitions(electrode_grid_shifted(:), 2);
z1 = electrode_channel_poitions(electrode_grid_unshifted(:), 3);
z2 = electrode_channel_poitions(electrode_grid_shifted(:), 3);

dists = sqrt((x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2);
mean_electrode_dist = mean(dists);


scouts_dist_per_electrode_dist = scouts_mean_distances / mean_electrode_dist;
scouts_dist_per_std = scouts_mean_distances ./ scouts_stds;


end

