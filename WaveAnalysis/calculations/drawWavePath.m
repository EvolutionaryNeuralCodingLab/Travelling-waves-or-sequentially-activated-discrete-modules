function waveCenterPath = drawWavePath(singleCrossings,singleHilbertAmps,startEndWave,En,varargin)
%DRAWWAVEPATH Draws the path of the wave by calculating in each frame the 
%"center of mass" of all crossings WITHIN WAVE, normalized by temporal 
%distance and hilbert amplitude. Will only look at first crossing within
%wave
%INPUT:
%   singleCrossings (channelsXcrossings)
%       All the crossings times. One of the four arrays that are recived in
%       'crossings' by the getHilbertCrossings function
%   singleHilbertAmps
%       Hilbert amplitude of the crossings. One of the four arrays that are 
%       recived in 'hilbertAmps' by the getHilbertCrossings function
%   startEndWave (1x2) 
%       Array with the indices which define the start and end of the 
%       window (reletive to beging of trial. i.e. 1 is the start of trial start).
%   En
%       Channel layout
%   Varargs (given as 'Name','Value' pairs):
%       flipEn: logical. Whether to flip En or not when when getting 
%               channel's position. Defult is 1 since 
%               convertChannelsToMovie flips En when converting recorded 
%               data (to match the En flipping of IntensityPhysicalSpacePlot).
%               If you're drawing for simulated data then 
%               convertChannelsToMovie is not used and therefor there
%               should be no flipping.
%       temporalWeightWidth: How wide is the temporal gaussian Sigma.
%           The maximal guassian sigma will be temporalWeightWidth times the
%           average distance between crossings (calculated from crossings
%           within startEndWindow). This weight is for calculating the COM
%           samples in the middle. Gaussian temporal weight decreases as 
%           you get closer to the edges. Default is 10 (Highley organized 
%           spatialy this is good, otherwise widen the width, e.g. change 
%           temporalWeightWidth to ~100)
%       normCoordinates - logical. Spatial coordinates will be normalized 
%           to the range [0,1] using maximal crossing position
%   
%OUTPUT:
%   waveCenterPath (nSamplesX2) 
%       coordinates (x,y) of the center of the wave.

flipEn=1;
temporalWeightWidth=10;
normCoordinates=1;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if flipEn
    En=flipud(En);
end

nSamples=startEndWave(2)-startEndWave(1)+1; %include edges samples in wave
nChannels=size(singleCrossings,1);

% %create arrays for all crossing's Amp, Position and times

crossAmps=[];
crossPos=[]; %(x,y)
crossTimes=[];

for i=1:nChannels
    [channelPosY,channelPosX]=find(En==i);
    crossInd=find(singleCrossings(i,:)>=startEndWave(1) & singleCrossings(i,:)<=startEndWave(2),1);
    if ~isempty(crossInd)
        crossTimes(length(crossTimes)+1)=singleCrossings(i,crossInd)-startEndWave(1);
        crossAmps(length(crossAmps)+1)=singleHilbertAmps(i,crossInd);
        crossPos(size(crossPos,1)+1,1:2)=[channelPosX,channelPosY];
    end
end

%create shorter temporal weighting width for boundary adjecent samples
meanCrossDiff=mean(diff(sort(crossTimes)));
temporalGaussWeightSigmas=[linspace(floor(temporalWeightWidth/5)*meanCrossDiff,temporalWeightWidth*meanCrossDiff,floor(nSamples/2)) linspace(temporalWeightWidth*meanCrossDiff,floor(temporalWeightWidth/5)*meanCrossDiff,ceil(nSamples/2))];

%calculate center of mass for each frame

waveCenterPath=zeros(nSamples,2);
for i=1:nSamples
   timeDiffs=abs(i-crossTimes); 
   tempWeight=exp(-timeDiffs.^2/(2*(temporalGaussWeightSigmas(i))^2));
   waveCenterPath(i,1)=sum(crossPos(:,1)'.*crossAmps.*tempWeight)/sum(crossAmps.*tempWeight);
   waveCenterPath(i,2)=sum(crossPos(:,2)'.*crossAmps.*tempWeight)/sum(crossAmps.*tempWeight);
end

%if normCoordinates, normalize coordinates reletive to chPos max
if normCoordinates
    waveCenterPath=waveCenterPath/max(crossPos(:));
end


end

