function [LFPWaveChannels,LFPWaveCrossingTimes] = getClusterFirstCrossings(clusterChannels,clusterTimes)
%GETCROSSINGSFROMCLUSTER recieves an array of all the channels participated
%in a crossing cluster clusterChannels (as recieved from getTrialClusters) 
%and their crossings times clusterTimes, and returnes just the first
%crossing per channel. LFPWaveChannels,LFPWaveCrossingTimes are ordered
%according to channel num

%%delete crossings other than the first for each channel
[LFPWaveCrossingTimes,orderOfCrossings]=sort(clusterTimes); %This is to make sure that the first crossings of a channel will be chosen
LFPWaveChannels=clusterChannels(orderOfCrossings);
[~, ind] = unique(LFPWaveChannels);
repeat_ind = setdiff(1:length(LFPWaveChannels), ind);
LFPWaveChannels(repeat_ind)=[];
LFPWaveCrossingTimes(repeat_ind)=[];
%reorder according to channel nums
[LFPWaveChannels,orderOfChannels]=sort(LFPWaveChannels);
LFPWaveCrossingTimes=LFPWaveCrossingTimes(orderOfChannels);

    
end

