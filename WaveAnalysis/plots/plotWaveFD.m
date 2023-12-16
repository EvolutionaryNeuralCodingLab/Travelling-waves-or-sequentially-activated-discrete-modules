function [h] = plotWaveFD(time,FD,waveChannels,waveTimes,pathChannels,samplingFrequency,varargin)
%PLOTWAVEFD Plots the filtered data FD of examplary channels through which
%the waveCenterPath goes through. 
%   INPUT:
%   - time: ms. Time of FD from start of the trial
%   - FD (nChannelsX1XnSamples): Filtered data. Make sure contains only one
%   trial!
%   - waveCenterPath (nWaveSamplesX2): Coordinates for the wave center
%   path, from the first crossing to the last (this is given by
%   drawWavePath function)
%   - En - electrode layout
%   - waveChannels,waveTimes - channels and crossings times of wave
%   pattern, as given by getClusterFirstCrossings or getRelevantWaveTimes

scatterSize=75;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if size(FD,2)>1
   warning('FD array contains more than one trial!') 
end

f=figure;
nRepresentativeChannels=size(pathChannels,1);
representativeChannelsFD=squeeze(flipud(FD(pathChannels,1,:)))';
baseCumSum=[0 cumsum(max(representativeChannelsFD(:,1:end-1))*1.5)];
lines=plot(time,representativeChannelsFD+baseCumSum,'LineWidth',1.5);
cm=jet(nRepresentativeChannels);
for i=1:nRepresentativeChannels
    set(lines(i),'color',cm(nRepresentativeChannels-i+1,:));
end

hold on
for i=1:nRepresentativeChannels
    crossingSample=waveTimes(waveChannels==pathChannels(nRepresentativeChannels+1-i));
    crossingHeight=baseCumSum(i)+representativeChannelsFD(crossingSample,i);
    scatter(crossingSample/samplingFrequency*1000,crossingHeight,scatterSize,cm(nRepresentativeChannels+1-i,:),'filled')
end

%     pos=get(gcf,'Position')
% set(gcf,'Position',[681   635   380   380]);
axis square
xlabel('Time [ms]','FontSize',20)
ylabel('Filtered Data','FontSize',20)
set(gca,'ytick',[])
set(gca,'yticklabel',[])


if nargout>0
   h=f; 
end

end

