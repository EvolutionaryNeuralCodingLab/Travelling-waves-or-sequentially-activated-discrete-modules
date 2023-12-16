function movieData = convertChannelsToMovie(chData,En,varargin)
%convertChannelsToMovie converts a ChxSamples matrix to 
%HeightXWidthXnFrames movie matrix according to En channel map.
%NOTICE that En is flipped upside down. If for example En=[1 2;3 4]
%then data from channel 1 will appear at the bottm ledt corner.
% varargins:
%   BGVal (1x1 double)
%       Value of pixels with no channels assigned to them (Usually
%       corners). Default value is min(chData(:)) (or min of Channels)
%   Channels (1XnCh):
%       Channels to export - all other channels will contain BGVal
%   flipEn (logical) - if false En is not flipped. Defualt true.

flipEn=1;

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end


%"remove" bad channels and set background
if exist('Channels','var')
    badChannels=setdiff(En(:),Channels);
    badChannels=badChannels(~isnan(badChannels));
    if exist('BGVal','var')
       chData(badChannels,:)=BGVal;
    else
       BGVal=min(min(chData(Channels,:)));
       chData(badChannels,:)=BGVal;
    end
else
    if ~exist('BGVal','var')
       BGVal=min(chData(:));
    end
end



nCh=max(En(:));

[frameHeight,frameWidth]=size(En);

movieData=BGVal*ones(frameHeight,frameWidth,size(chData,2));

if flipEn
    En=flipud(En);
end
for i=1:nCh
    [chPosY,chPosX]=find(En==i);
    if size(chPosY)>0
        movieData(chPosY,chPosX,:)=chData(i,:);
    end
end

end

