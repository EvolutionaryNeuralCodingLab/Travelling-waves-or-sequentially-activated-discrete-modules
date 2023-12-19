function [ch_nums_ordered] = order_channels_by_crossing_times(crossings,start_end_wave,ch_nums)
%ORDER_CHANNELS_BY_CROSSING_TIMES finds the first phase crossing within
%start_end_wave and returns the order of the channels accordnig to it. If
%start_end_wave is not given it just finds the first phase.
% Temporal units are assumed to be in samples


if ~exist('start_end_wave','var')
    start_end_wave = [1 max(crossings(:))];
end

if ~exist('ch_nums','var')
    ch_nums = 1:size(crossings,1);
end

temp_crossings = crossings;
%zero out everything that's outside of range
temp_crossings((temp_crossings<start_end_wave(1) & temp_crossings>0) | temp_crossings>start_end_wave(2)) = 0;
%find first non-zero
[has_max_in_wave, max_index] = max( temp_crossings ~=0, [], 2 );
%drop bad channels
max_index = max_index(has_max_in_wave);
ch_nums = ch_nums(has_max_in_wave);
temp_crossings = temp_crossings(has_max_in_wave,:);
%get max per raw; using 1:size(temp_crossings,1) and not ch_nums on purpose
first_crossings_times = temp_crossings(sub2ind(size(temp_crossings),1:size(temp_crossings,1),max_index'));
% now sort
[~, I] = sort(-first_crossings_times);

ch_nums_ordered = ch_nums(I);

end

