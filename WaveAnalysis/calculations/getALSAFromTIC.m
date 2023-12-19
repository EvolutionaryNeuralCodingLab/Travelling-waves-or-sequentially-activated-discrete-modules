function ALSA = getALSAFromTIC(ticPath,startTime,window_ms,En,samplingFrequency,varargin)
%GETALSAPEAKSFROMTIC Calculates ALSA (average local spiking activity) from
%startTime to startTime+window_ms. If length(startTime)>1, getALSAFromTIC
%returns the average across trials
%   INPUT:
%       - ticPath: full path to .mat file containing recording's t,ic 
%       - startTime (ms) (1XnTrials) - Time from the recording begining for 
%       which to calculate ALSA. If an array, getALSAFromTIC returns
%       average over trials
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
%           - nNearestChannels - number of nearest neighbors to calc ALSA
%           over (4 or 8). Default is 4.
% 
%   OUPUT
%       - ALSA (nChannelsXnSamples) - the average ALSA over trials.
%       nChannels taken from En. nSamples calculated according to startEndWave

spikeRateSmoothing=2000;
startEndWave=[1 window_ms*samplingFrequency/1000];
nNearestChannels=4;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

window_samples=startEndWave(2)-startEndWave(1)+1;
nTrials=length(startTime);
nCh=max(En(:));

ALSA=zeros(nCh,window_samples);

for i=1:nTrials

    %calc and smooth spiking rate
    spikingRateDataFormat=tic2FireRate(ticPath,startTime(i),window_ms,En,samplingFrequency,'outputFormat','dataFormat','slidingWindowSize',spikeRateSmoothing);
    spikingRateDataFormatSmoothed=smoothdata(spikingRateDataFormat','gaussian',spikeRateSmoothing)';

    %calculate ALSA
    ALSA=ALSA+calcALSA(spikingRateDataFormatSmoothed,'En',En,'nNearestChannels',nNearestChannels);
end

ALSA=ALSA/nTrials;

end

