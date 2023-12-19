function [f] = plotSingleHilbertCrossing(singleCrossings,crossingsAmps,FD,crossingType,channelShown,varargin)
%PLOTSINGLEHILBERTCROSSING plots the occurences of a specific crossings,
%with color map indicating the hilbert amplitude. It also plots filtered
%data from a single exemplary channel
%   singleCrossings are the specific crossings (upwards, downwards,
%   inhibition or excitation) as recived by getHilbertCrossings.
%   crossingsAmps are the amplitudes at these times. 
%   FD (1XnSamples) is the exemplary filtered data from a specific channel
%   crossingType (string) is the crossing type (minima/inhibition etc.)
%   channelShown is the channel number of the data being shown (FD)
%   Varargins (given as 'Key'Value pairs:
%       - Spikes (logical nChXnSamples)
%      Plots a red circle where Spikes is true
%       - spikesStyle (String)
%      Style of spikes. Default is '.r'
%       - Title (string)
%      Figure title
%       - clusterLimits (nClusterX2)
%      Plot also crossings clusters
%       - plotLegend (logical)
%      True to plot legend (default), False for no legend.
%      moreChannelsData overrides this to be false.
%       - plotColorbar: default is true
%       - CrossingsVerticalOffset (1x1)
%      Move Crossings scatter up or down. Default is 0.
%       - moreChannelsData (nChXnSamples)
%      Plot additional data in figure. Each channel will be normalized
%      and plotted at the hight of the row number (channel i will be 
%      below channel i+1). If this is given, plotLegend is set to
%      false. Whoever's reading this is welcome to fix the legend it
%      creates - I didn't have time.
%       - moreChannelsNumbers (nCX1) configures the vertical position of
%      every channel in moreChannelsData (allows to plot less channels then
%      total channel numbers, or to order it in different way).
%       - normalizeChannels (1x1 logical) Normelaize each channel by twice
%      its mean. default is true.
%       - timeInms - logical. Display crossings times in ms, using
%       sample2ms. Default is 0.
%       - sample2ms - 1000/samplingFrequency
%       - sz - crossings' circle marker size. Default 25
%       - spikeMarkerSize - default 6
%       - minHilbertAmp - only show crossings with whose hilbert amplitude
%       is above minHilbertAmp. Default is 0
%       - order - order in which to plot channels

%%QA make sure that ms2sample and timeInms works

plotColorbar=true;
plotLegend=true;
CrossingsVerticalOffset=0;
plotAdditionalChannels=0;
normalizeChannels=1;
spikesStyle='.r';
plotSpikes=0;
timeInms=0; 
spikeMarkerSize=6;
sz=25;
minHilbertAmp=0;
order=1:size(singleCrossings,1);

f=figure;
hold on

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if exist('Spikes','var')
    plotSpikes=1;
end

if exist('moreChannelsData','var')
    plotAdditionalChannels=1;
    plotLegend=0;
    if ~exist('moreChannelsNumbers','var')
        moreChannelsNumbers=(1:size(moreChannelsData,1))';
    end
end

chNum=1:size(singleCrossings,1);

[~,chOrder]=sort(order);

if plotSpikes
    for i=chNum
       spikeInd=find(Spikes(i,:));
       if timeInms
           if ~exist('sample2ms','var')
              error('To plot time in ms instead of seconds sample2ms must be given') 
           end
       end
       plot(spikeInd*sample2ms, chOrder(i)*ones(1,length(spikeInd)),'.k','markerSize',spikeMarkerSize);
    end
end


%get rid of low hilbert amps
singleCrossings(crossingsAmps<=minHilbertAmp)=0;
crossingsAmps(crossingsAmps<=minHilbertAmp)=0;

h(1)=scatter(singleCrossings(chNum(1),singleCrossings(chNum(1),:)~=0)*sample2ms,chOrder(chNum(1))*ones(1,numel(singleCrossings(chNum(1),singleCrossings(chNum(1),:)~=0)))+CrossingsVerticalOffset,sz,squeeze(crossingsAmps(1,singleCrossings(chNum(1),:)~=0))');
hold on
for i=chNum(2:end)
        scatter(singleCrossings(i,singleCrossings(i,:)~=0)*sample2ms,chOrder(i)*ones(1,numel(singleCrossings(i,singleCrossings(i,:)~=0)))+CrossingsVerticalOffset,sz,crossingsAmps(i,singleCrossings(i,:)~=0)'); %make crossingsAmp column vec to avoid RGB syntex in case there is only 3 element
end
% h(2)=plot(1:length(FD)*sample2ms,'b');
% h(3)=plot((0:numel(FD))*sample2ms,channelShown*ones(1,numel(FD)+1)+CrossingsVerticalOffset,'--k');

h(2)=plot((1:length(FD))*sample2ms,FD,'b');
if ~isempty(channelShown)
    h(3)=plot(0:numel(FD),channelShown*ones(1,numel(FD)+1)+CrossingsVerticalOffset,'--k');
end


if plotLegend
%legend hack
h(1)=plot(nan,nan,'ob');
h(2)=plot(nan,nan,'b');
h(3)=plot(nan,nan,'--k');
h(4)=plot(nan,nan,'.r');
end

if plotSpikes
    if plotLegend
        legend([h(1) h(2) h(3) h(4)],{crossingType,'Filtered Data',['Current Channel (' num2str(channelShown) ')'],'Spikes'})
    end
else
    if plotLegend
        legend([h(1) h(2) h(3)],{crossingType,'Filtered Data','Current Channel'})
    end
end


if exist('Title','var')
    title(Title)
end

if exist('clusterLimits','var')
    nClusters=size(clusterLimits,1);
    for i=1:nClusters
        plot(clusterLimits(i,:)*sample2ms,[0 0],'LineWidth',2,'Color','k')
        if plotLegend
        %remove added legend
            hLegend = findobj(gcf, 'Type', 'Legend');
            newLegend=hLegend.String(1:end-1);
            legend(newLegend)
        end
    end
end

if plotColorbar
    hcb=colorbar;
    title(hcb,'Hilbert Amplitude [uV]');
end
if timeInms
    xlabel('Time [ms]')
else
    xlabel('Time [Samples]')
end

% ylabel(['Filtered Data (Channel' num2str(channelShown) ') [uV]'])
ylabel('Channels')
    

if plotAdditionalChannels
    if size(moreChannelsNumbers,2)>1
        warning('moreChannelsNumbers size should be nChX1!')
    end
    moreDataSamples=size(moreChannelsData,2);
    %first normalize
    if normalizeChannels
        moreChannelsData=moreChannelsData./(repmat(std(moreChannelsData,0,2),1,moreDataSamples));
    end
    %add relevant height
    moreChannelsData=moreChannelsData+repmat(moreChannelsNumbers,1,moreDataSamples);
%     hLegend = findobj(gcf, 'Type', 'Legend');
%     keepLegend=hLegend.String;
%     nLegends=numel(hLegend.String);
    plot((1:size(moreChannelsData,1))*sample2ms,moreChannelsData');
%     hLegend = findobj(gcf, 'Type', 'Legend');
%     keepLegend=hLegend.String;
%     legend(keepLegend(1:nLegends))
end

end


