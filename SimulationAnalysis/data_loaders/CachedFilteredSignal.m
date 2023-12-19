classdef CachedFilteredSignal
    %CACHEDFILTEREDSIGNAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cached_folder
        cached_prefix
    end
    
    properties (Constant)
        DEFAULT_CACHED_FOLDER = [get_analysis_code_base_path 'cached_filtered_data/'];
        DEFAULT_CACHED_PREFIX = 'electrode_signal_bandpass_filtered_';
    end

    methods
        function obj = CachedFilteredSignal(cached_folder,cached_prefix)
            %CACHEDFILTEREDSIGNAL Construct an instance of this class
            %   Detailed explanation goes here
            if exist('cached_folder','var')
                obj.cached_folder = cached_folder;
            else
                obj.cached_folder = CachedFilteredSignal.DEFAULT_CACHED_FOLDER;
            end
            if exist('cached_prefix','var')
                obj.cached_prefix = cached_prefix;
            else
                obj.cached_prefix = CachedFilteredSignal.DEFAULT_CACHED_PREFIX;
            end

        end
        

        function hached_full_path = get_hached_full_path(obj, data)
            H = obj.get_data_hach(data);
            hached_full_path = [obj.cached_folder obj.cached_prefix H '.mat'];
        end
    end
    
    methods(Static)
        function filtered_signal_cache_struct = create_filtered_signal_cache_struct(electrode_signals,ws, fs, crop_at, pad_samples, flip_signal,n_reps)
            if nargin==3 
                filtered_signal_cache_struct.electrode_signals = electrode_signals;
                filtered_signal_cache_struct.ws = ws;
                filtered_signal_cache_struct.fs = fs;
            elseif nargin==7 
                filtered_signal_cache_struct.crop_at = crop_at;
                filtered_signal_cache_struct.pad_samples = pad_samples;
                filtered_signal_cache_struct.flip_signal = flip_signal;
                filtered_signal_cache_struct.n_reps = n_reps;
            else 
                error('wrong number of input parameters')
            end
        end
    end
    methods (Hidden)
        function H = get_data_hach(obj, data)
            Engine = java.security.MessageDigest.getInstance('MD5');
            H = obj.core_hach(data, Engine);
            H = sprintf('%.2x', H);   % To hex string
        end

        function H = core_hach(obj, data, engine)
            % Consider the type of empty arrays:
            S = [class(data), sprintf('%d ', size(data))];
            engine.update(typecast(uint16(S(:)), 'uint8'));
            H = double(typecast(engine.digest, 'uint8'));
            if isa(data, 'struct')
               n = numel(data);
               if n == 1  % Scalar struct:
                  F = sort(fieldnames(data));  % ignore order of fields
                  for iField = 1:length(F)
                     H = bitxor(H, obj.core_hach(data.(F{iField}), engine));
                  end
               else  % Struct array:
                  for iS = 1:n
                     H = bitxor(H, obj.core_hach(data(iS), engine));
                  end
               end
            elseif isempty(data)
               % No further actions needed
            elseif isnumeric(data)
               engine.update(typecast(data(:), 'uint8'));
               H = bitxor(H, double(typecast(engine.digest, 'uint8')));
            elseif ischar(data)  % Silly TYPECAST cannot handle CHAR
               engine.update(typecast(uint16(data(:)), 'uint8'));
               H = bitxor(H, double(typecast(engine.digest, 'uint8')));
            elseif iscell(data)
               for iS = 1:numel(data)
                  H = bitxor(H, obj.core_hach(data{iS}, engine));
               end
            elseif islogical(data)
               engine.update(typecast(uint8(data(:)), 'uint8'));
               H = bitxor(H, double(typecast(engine.digest, 'uint8')));
            elseif isa(data, 'function_handle')
                H = bitxor(H, obj.core_hach(functions(data), engine));
            else
               warning(['Type of variable not considered: ', class(data)]);
            end
        end
    end
end

