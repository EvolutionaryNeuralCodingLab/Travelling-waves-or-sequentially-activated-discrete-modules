function [relevantChannels,relevantCrossingTimes,relevantALSATimes,irelevantChannels] = getRelevantWaveTimes(LFPchannels,LFPtimes,ALSAchannels,ALSA_Locs,nCh)
%GETRELEVANTTIMES Summary of this function goes here
%   Detailed explanation goes here

[LFPWaveChannels,LFPWaveCrossingTimes]=getClusterFirstCrossings(LFPchannels,LFPtimes);

[relevantChannels,LFPind,ALSAind]=intersect(LFPWaveChannels,ALSAchannels);

relevantCrossingTimes=LFPWaveCrossingTimes(LFPind);
relevantALSATimes=ALSA_Locs(ALSAind);

irelevantChannels=setdiff(1:nCh,relevantChannels);

    
end

