function [PLM,f,ax] = plotCrossingsPhysical(selectedCrossings,startEndWave,En,hilbertAmps,varargin)
%plotCrossingsPhysical uses IntensityPhysicalSpacePlot to plot the
%time to the first crossing in selectedCrossings for each channel. It looks
%at the window defined by startEndWave. t=0 is defined by the time that the 
%first channel crossing within the window has a crossing.
%Values of selectedCrossings and startEndWave must have the same units and
%starting time.
%INPUT:
%   selectedCrossings (channelsXcrossings)
%       All the crossings times
%   startEndWave (1x2) 
%       Array with the indices which define the start and end of the 
%       window in which the first crossing is looked for.
%   En
%       Channel layout
%   Possible Varargs (given as 'key',value pairs):
%   Units (string)
%       Units of crossing times and startEndWave. Default is ms.
%   Title (string)
%       Figure title
%   plotElectrodeNumbers (logical)
%      Default is false
%   plotColorBar (logical)
%       Default is true
%OUTPUT:
%   PLM (1XnCh)
%       Phase Latency Map - the times all channel reached the crossings
%       for the first time, in same units as selectedCrossings. The time of
%       first channel to cross is defined as t=0. Channels which did not 
%       cross will contain NaNs.
%TODO: Add option not to send hilbertAmps

f=figure;
Units='ms';
plotElectrodeNumbers=false;
plotSizeBar=0;
plotColorBar=1;
plotGridLines=0;
figSizeUnits='centimeters';
axesPosition=[2.14,0.34,0.19,1.9];
hCbarFontSize=12;
% hCbarPosition=[13.8113    8.4667    5.9267    2.6458];
hCbarPosition=[12.3826,7.8052,0.7143,2.4871];
hCbarYlabel={'Latency','[ms]'};
% hCbarTextPosition=[-4,0.5,0];
hCbarTextPosition=[1.26,0.52,0];
hCbarFontUnits='points';
hCbarYlablFontSize=8;
markerSize=18;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

chNum=1:size(selectedCrossings,1);

pT=[];
hAmp=[];
channels=[];
for i=chNum
    findCross=find(selectedCrossings(i,:)>=startEndWave(1) & selectedCrossings(i,:)<=startEndWave(2));
    if numel(findCross)==0
        pT(i)=NaN;
        if ~isempty(hilbertAmps)
            hAmp(i)=0;
        end
%         i
    else
%         if numel(findCross)>1 , i ,end %notify when a channel has more than 1 crossing within the window
        pT(i)=selectedCrossings(i,findCross(1))-startEndWave(1);
        if ~isempty(hilbertAmps)
            hAmp(i)=hilbertAmps(i,findCross(1));
        end
        channels(length(channels)+1)=i;
    end
end
firstCrossing=min(pT);

pT=(pT-firstCrossing);

%plot the times of crossing in physical space
if ~isempty(hilbertAmps)
    [hCbar,ax]=IntensityPhysicalSpacePlot(chNum,pT,En,'plotElectrodeNumbers',plotElectrodeNumbers,'markerSize',hAmp*60/max(hAmp),'plotSizeBar',plotSizeBar,'plotColorBar',plotColorBar,'plotGridLines',plotGridLines);
else
    [hCbar,ax]=IntensityPhysicalSpacePlot(chNum,pT,En,'plotElectrodeNumbers',plotElectrodeNumbers,'markerSize',markerSize,'plotColorBar',plotColorBar,'plotGridLines',plotGridLines);
end

% set(gcf,'Units',figSizeUnits,'Position',figPosition);
% set(gca,'Units',figSizeUnits,'Position',axesPosition);

if plotColorBar
    hCbar.Units=figSizeUnits;
    hCbar.FontSize=hCbarFontSize;
    hCbar.Position=hCbarPosition;
    ylabel(hCbar,hCbarYlabel,'Position',hCbarTextPosition,'FontUnits',hCbarFontUnits,'FontSize',hCbarYlablFontSize);
end
if exist('Title','var')
    title(Title)
end 

% xlabel(['Note: Channels with no crossings in window is assigned the value of ' num2str((nullValue-firstCrossing)*sample2ms) 'ms'],'Color',[1 1 1]*0.5,'FontSize',9)


if nargout>0
    PLM=pT;
end

end

