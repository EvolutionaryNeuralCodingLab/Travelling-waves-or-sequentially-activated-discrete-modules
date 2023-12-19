function crop_at = get_crop_at_if_empty(crop_at,signal_length)
if isempty(crop_at)
    crop_at = [1, signal_length];
end

end

