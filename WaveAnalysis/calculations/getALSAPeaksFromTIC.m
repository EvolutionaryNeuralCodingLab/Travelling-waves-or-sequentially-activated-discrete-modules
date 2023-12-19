function [ALSA_Locs,relevantChannels] = getALSAPeaksFromTIC(ticPath,startTime,window_ms,En,samplingFrequency,varargin)
%GETALSAPEAKSFROMTIC Calculates ALSA (average local spiking activity) from
%startTime to startTime+window_ms and returns location of each channel's 
%first peak.
%   INPUT:
%       - ticPath: full path to .mat file containing recording's t,ic 
%       - startTime (ms) - Time from the recording begining for which to
%       calculate ALSA.
%       - window_ms (ms) - the window length
%       - En - electrode layout
%       - samplingFrequency - sample/s
%       - Possible Varargins: (given as 'Key','Value' pairs)
%           - spikeRateSmoothing (samples) - The window by which the 
%           averaging and smoothing is done. Default is 2000 (100ms for 
%           20k samling rate)
%           - startEndWave (samples) - Subwindow from which to get first 
%           ALSA maxima. i.e. First all alsa maxima are found, and then the
%           first maximum within this window is taken. Default is 
%           [1 nSamples]
%           - outputUnits - 'samples' (default) or 'ms'
%           - onsetType - Method for detecting local spiking onset. Can be
%           either 'firstMax' (default) or 'stdCrossings'.
%           - nSTD - relvent for onsetType='stdCrossings'. ALSA
%           onset will be detected when spiking rate crosses nSTD stds.
%           Default is 2
%           - nNearestChannels - number of nearest neighbors to calc ALSA
%           over (4 or 8). Default is 4.
%
%   OUPUT
%       - ALSA_Locs (1XnRelevantChannels) - locations of first ALSA maxima, given in samples
%       (couting from the startTime, not startEndWindow(1)).
%       - relevantChannels (1XnRelevantChannels) - List of the channels 
%       that had an ALSA max within window

onsetType='firstMax';
nSTD=2;
outputUnits='samples';
spikeRateSmoothing=2000;
startEndWave=[1 window_ms*samplingFrequency/1000];
nCh=max(En(:));
nNearestChannels=4;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

%calc ALSA
ALSA = getALSAFromTIC(ticPath,startTime,window_ms,En,samplingFrequency,'spikeRateSmoothing',spikeRateSmoothing,'nNearestChannels',nNearestChannels);

% Find each channel's first peak withihn window
ALSA_Locs=zeros(1,nCh);

if strcmp(onsetType,'firstMax')
    for i=1:nCh
        [~,chlocs]=findpeaks(ALSA(i,:));
    %     firstPeakInWindow=find(chlocs>=startEndWave(1) & chlocs<=startEndWave(2),1);
        if isempty(chlocs) || chlocs(1)<startEndWave(1) || chlocs(1)>startEndWave(2)
            ALSA_Locs(i)=NaN;
        else
            ALSA_Locs(i)=chlocs(1);
        end
    end
else 
    chSTDs=std(ALSA,0,2);
    for i=1:nCh
        chLoc=find(ALSA(i,1:end-1)<nSTD*chSTDs(i)&ALSA(i,2:end)>=nSTD*chSTDs(i),1);
        if isempty(chLoc) || chLoc<startEndWave(1) || chLoc>startEndWave(2)
            ALSA_Locs(i)=NaN;
        else
            ALSA_Locs(i)=chLoc;
        end
    end
end

relevantChannels=find(~isnan(ALSA_Locs));
irrelevantChannels=setdiff(1:nCh,relevantChannels);
ALSA_Locs(irrelevantChannels)=[];

if strcmp(outputUnits,'ms')
   ALSA_Locs=ALSA_Locs/samplingFrequency*1000;
end
end

