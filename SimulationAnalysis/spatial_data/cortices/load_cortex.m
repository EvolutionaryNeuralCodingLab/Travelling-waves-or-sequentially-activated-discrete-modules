function [cortex] = load_cortex(cortex_name)
%LOAD_CORTEX_VERTICES returns a cortex struct given a cortex_name
%   Cortex struct contains two fields:
%       cortex.name - name of cortex (taken from cortex_name
%       cortex.verices_positions - 3d positions of cortex vertices. Will be
%                                  used mainly for scout positioning.

% %don't know why, dictionary didn't work
% cortex_name_to_filename = dictionary("306716V", ...
%                                      "cortex306716V.mat",...
%                                      "15002V",...
%                                      "cortex_15002V.mat.mat"...
%                             );
% Going with cell array for now

cortex_name_to_filename = {'306716V',...
                           'cortex306716V.mat',...
                           '15002V',...
                           'cortex_15002V.mat',...
    };

ind = find(contains(cortex_name_to_filename,cortex_name));

if isempty(ind)
    error(strcat("Didn't find the cortex you were looking for ('", cortex_name, "'). Available Cortices Names: ", strjoin(cortex_name_to_filename(1:2:end), ', ')))
elseif length(ind)<=2
    cortex_vars = load(cortex_name_to_filename{ind(1)+1});
    field_name = fieldnames(cortex_vars);
    field_name = field_name{1};
    cortex.vertices_positions = cortex_vars.(field_name).Vertices;
else
    error(strcat("Cortex name ('", cortex_name, "') appeas more than once. Please be more specific. Requested name found in: ", strjoin(cortex_name_to_filename(ind), ', ')))
end
cortex.name = cortex_name;

end

