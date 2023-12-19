function filled_plotting_config = fill_missing_configs_with_defaults(plotting_config, default_plotting_config)
%FILL_MISSING_CONFIGS_WITH_DEFAULTS fill plotting_config with defaults
% defined in default_plotting_config that are missing from plotting_config.
% If default_plotting_config is not given, get_default_plotting_config() is
% used.
%   plotting_config controls the way analyze_simulation_from_table is used.
%   Along with what plots to create, it controls filtering, noise addition,
%   etc.

if ~exist('default_plotting_config','var')
    default_plotting_config = get_default_plotting_config();
end
filled_plotting_config = plotting_config;
config_fields = fieldnames(filled_plotting_config);
all_fields = fieldnames(default_plotting_config);
for i=1:numel(all_fields)
    if ~any(strcmp(config_fields,all_fields{i}))
        filled_plotting_config.(all_fields{i}) = default_plotting_config.(all_fields{i});
    end
end

