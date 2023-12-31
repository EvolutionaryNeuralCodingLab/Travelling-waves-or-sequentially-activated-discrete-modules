function [electrode_channel_poitions, cortex, simulation_table, electrode_grid, save_path] = apply_path_configurations_and_load(path_to_channel_positions,cortex_name,path_to_results_table,electrode_grid_name, save_path)
%APPLY_PATH_CONFIGURATIONS Sets and loads parameters use for analysis, 
% according to chosen paths and configurations

%load channel positions
electrode_channel_poitions = load(path_to_channel_positions );
electrode_channel_poitions = electrode_channel_poitions.channel_locations;

%load cortex struct
cortex = load_cortex(cortex_name);

%load simulation table
simulation_table = load(path_to_results_table);

%load electrode grid (En)
base_path = get_analysis_code_base_path();
electrode_grids_path = [base_path 'spatial_data/elecrode_grids/'];
electrode_grid_path = [electrode_grids_path,electrode_grid_name,'.txt']; 
if ~isfile(electrode_grid_path)
    folder_files = dir(electrode_grids_path);
    folder_filenames = extractfield(folder_files,'name');
    available_grids_names = cellfun(@(x) x(1:end-4),folder_filenames(3:end),'UniformOutput',false);
    error(strcat("Selected electrode grid ('", electrode_grid_name, "') Does not exists. Please choose one of the following: '", strjoin(available_grids_names, "', '"),"'"))
end
electrode_grid = readmatrix(electrode_grid_path);

%apply electrode grid name to save path
if ~isempty(save_path)
    save_path = [save_path, electrode_grid_name, '/'];
    if ~exist(save_path, 'dir')
       mkdir(save_path)
    end
end
end