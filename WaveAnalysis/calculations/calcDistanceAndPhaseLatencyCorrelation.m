function [r,distances] = calcDistanceAndPhaseLatencyCorrelation(channels,crossingTimes,En)
%CALCDISTANCEANDPHASELATENCYCORRELATION Calculates the correlation between
%each channel's phase crossing latency and the euclidian distance between 
%the channel and the first channel that reached the phase
%
%   INPUT:
%       - channels - list of channels (like the ones recived by
%       getTrialClusters)
%       - crossingTimes - phase latencies of 'channels' (like the ones 
%       recived by getTrialClusters)
%       - En - electrode layout
%	OUTPUT:
%       - r - correlation coefficient between distance and crossing time
%       - distances - distances between the channels' position to the
%       position of the first channel that crossed. distances(i)
%       corresponds to the channel in crossingTimes

chPos = calcChannelsPosition(En);
%leave only relevant channels

%get first channel to cross
waveStartCh=channels(find(crossingTimes==min(crossingTimes),1));
%calc dists from the channel

dists=sqrt((chPos(:,1)-chPos(waveStartCh,1)).^2+(chPos(:,2)-chPos(waveStartCh,2)).^2);

distances=dists(channels);

r=corrcoef(distances,crossingTimes);

r=r(1,2);

end

