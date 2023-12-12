function [path] = get_analysis_code_base_path()
%GET_ANALYSIS_CODE_BASE_PATH returns a string with the simulation analysis 
% code base

parts = strsplit(fileparts(mfilename('fullpath')), filesep);
path = [strjoin(parts(1:end-1), filesep) filesep];


end

