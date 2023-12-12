function crop_at = get_crop_at_if_empty(crop_at,signal_length)
%GET_CROP_START_END_IF_EMPTY Summary of this function goes here
%   Detailed explanation goes here
if isempty(crop_at)
    crop_at = [1, signal_length];
end

end

