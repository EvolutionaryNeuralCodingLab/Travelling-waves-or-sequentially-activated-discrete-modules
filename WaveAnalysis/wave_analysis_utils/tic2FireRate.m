function spikingRate = tic2FireRate(ticPath,startTime,window_ms,En,samplingFrequency,varargin)
%tic2FireRate calculates the spiking rate (spike/s) for each channel
%starting from startTime to startTime+window_ms and returns it as 3d 
%matrix: (frameHeightXFrameWidthXFrames) according to channel layout in En
%If length(startTime)>1, tic2FireRate returns the average across trials
%   Input:
%       -   ticPath: The path to t,ic containg saved mat
%       -   startTime (ms): 1XnTrials array with the start time of the 
%           wave. If nTrials>1 spikingRate returns the average over trials
%       -   En: Electrode layout
%       -   Varargs (given as 'Name','Value' pairs):
%       	-   slidingWindowSize (1X1): size of the moving window in units
%       	of samples. Default is 10000 (0.5s for 20kHz)
%           -   outputFormat (string): 
%               - "movieFormat": frameHeightXFrameWidthXFrames. Default.
%               - "dataFormat": Output is nChXnSamples (like the output of 
%                   recordingObject.getData). 
%
%   Output:
%       - spikingRate: First two dimension are the same as En. Third
%       dimention is time in units of samples (with the length of the
%       wave). If outputFormat is set to dataFormat, dimensions are
%       nChXnSamples
%
%   TODO: - Add varargin so that slidingWindowSize will be in different units
%         - MAKE SURE DIMENSIONS ARE frameHeightXFrameWidthXFrames (or,
%         fix description).

slidingWindowSize=10000; %samples
outputFormat='movieFormat';

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

nCh=max(En(:));
nSamples=window_ms*samplingFrequency/1000;
nTrials=length(startTime);

spikingRate=zeros(nCh,nSamples);

for i=1:nTrials
    binSpikes = getSpikeBinMatByChannel(ticPath,startTime(i),startTime(i)+window_ms,samplingFrequency,max(En(:)));

    spikingRate = spikingRate+binSpikes2fireRate(binSpikes,samplingFrequency,'slidingWindowSize',slidingWindowSize);
end

spikingRate=spikingRate/nTrials;


if strcmp(outputFormat,'movieFormat')
    spikingRate=convertChannelsToMovie(spikingRate,En);
end
     
end

