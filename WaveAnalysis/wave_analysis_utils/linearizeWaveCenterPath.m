function [roundWaveCenterPath,centerMoveTimes,pathChannels,waveChannelsPos] = linearizeWaveCenterPath(waveCenterPath,En)
%LIEARIZEWAVECENTERPATH recieve the wave center path coordinates and
%calculates the linear path starting from waveCenterPath(1,:) to
%waveCenterPath(2,:). It returns just the rounded linear path, the unique 
%channels it goes through and the times it went through them
%   INPUT:
%   - waveCenterPath - nSamplesX2 (x,y coordinates as recieved by drawWavePath
%   - En - electrode layout
%   OUTPUT:
%   - roundWaveCenterPath - 
%   - centerMoveTimes - the samples in which the center has moved
%   - pathChannels - the unique channel numbers through which the center
%   has went through

nSamples=size(waveCenterPath,1);
linearPathFromStartToEnd=[linspace(waveCenterPath(1,1),waveCenterPath(end,1),nSamples)' linspace(waveCenterPath(1,2),waveCenterPath(end,2),nSamples)'];
roundWaveCenterPath=round(linearPathFromStartToEnd); %linear path
centerMoveTimes=find(diff(roundWaveCenterPath(:,1)')~=0 | diff(roundWaveCenterPath(:,2)')~=0);
waveChannelsPos=roundWaveCenterPath(centerMoveTimes+1,:);

pathChannels=En(sub2ind(size(En),size(En,1)+1-waveChannelsPos(:,2),waveChannelsPos(:,1))); %flipped ud



end

