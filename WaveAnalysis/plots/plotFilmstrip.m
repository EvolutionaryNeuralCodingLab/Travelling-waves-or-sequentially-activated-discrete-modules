function [f] = plotFilmstrip(data,flimstripStartEnd,En,waveChannelsPos,samplingFrequency,varargin)
%PLOTFILMSTRIP plots a filmstrip of the data in data, at fixed intervals
%inside filmestripStartEnd. It also draws the elecrode positions described
%in waveChannelPos
%   INPUT:
%   - data (nChX1XnSamples) - the data to draw. MAKE SURE CONTAIS ONLY ONE
%   TRIAL!
%   - filmstripStartEnd (1X2) - start and end samples for the
%   filmstrip. 12 frames equally seperated are taken from this.
%   - En - electrode layout by which the frames are set up from data
%   - waveChannelsPos (nElectrodesX2) - x,y positions of electrodes to
%   plot. flipped ud.
%   - possible varargin:
%       colorBarWidthCM - width of colorbar in cm. default 0.25
%       figPos - 'Position' value of the figure in cm
%       plotTitleTimes - logical. Default is 0
%       colorbarFontSize - default is 12
%       boxEdgeSize - each frame edge size in cm. default 0.75
%       boxPadding - padding between subplots. default is 0.25

colorBarWidthCM=0.15;
figPos=[19,14,7.3,2.3];
plotTitleTimes=0;
colorbarFontSize=12;
boxEdgeSize=0.75; %cm
boxPadding=0.25; %cm

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if size(data,2)>1
   warning('Data array contains more than one trial!') 
end

f=figure;



numOfElectrodes=size(waveChannelsPos,1);
representativeTimes=round(linspace(flimstripStartEnd(1),flimstripStartEnd(2),12)); %12 snippets
movieData = convertChannelsToMovie(squeeze(data(:,1,:)),En);
nRep=length(representativeTimes);

n = ceil(nRep/2); % How many boxes in a row
% ax = ones(2,n); %Handle so you can plot afterward


minMovieVal=min(min(min(movieData(:,:,representativeTimes))));
maxMovieVal=max(max(max(movieData(:,:,representativeTimes))));


for i = 1:n
    ax(1,i) = axes('Parent',f,'Units','centimeters','Position',[boxPadding+(boxPadding+boxEdgeSize)*(i-1) 2*boxPadding+boxEdgeSize boxEdgeSize boxEdgeSize],'Visible','off','XTick',[],'YTick',[]);
    imagesc(movieData(:,:,representativeTimes(i)),[minMovieVal maxMovieVal])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    if plotTitleTimes
        title(['T=' num2str(round(representativeTimes(i)/samplingFrequency*1000)) 'ms'])
    end
    ax(2,i) = axes('Parent',f,'Units','centimeters','Position',[boxPadding+(boxPadding+boxEdgeSize)*(i-1) boxPadding boxEdgeSize boxEdgeSize],'Visible','off','XTick',[],'YTick',[]);
    imagesc(movieData(:,:,representativeTimes(i+n)),[minMovieVal maxMovieVal])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])

    if plotTitleTimes
        title(['T=' num2str(round(representativeTimes(i)/samplingFrequency*1000)) 'ms'])
    end
end

% imagesc(ax(1,1),movieData(:,:,representativeTimes(1)),[minMovieVal maxMovieVal])
ax2 = axes;
P = get(ax(1,1),'Position');
XLIM = get(ax(1,1),'XLim');
YLIM = get(ax(1,1),'YLim');
PA = get(ax(1,1),'PlotBoxAspectRatio');
set(ax2,'Units','centimeters','Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
scatter(ax2,waveChannelsPos(:,1),size(En,1)+1-waveChannelsPos(:,2),10,linspace(0,1,numOfElectrodes),'filled')
%%Link them together
linkaxes([ax(1,1),ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

%%Give each one its own colormap
colormap(ax(1,1),'parula')
colormap(ax2,'jet')
% %%Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
% cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);


set(f,'Units','centimeters','Position',figPos);
% set(subplot(2,ceil(nRep/2),i),'Units','centimeters')
hp4 = get( ax(2,end),'Position');
colorbar(ax(1),'Units','centimeters','Position', [hp4(1)+hp4(3)+0.1  hp4(2)  colorBarWidthCM  hp4(2)+hp4(3)*2.1],'FontSize',colorbarFontSize)



end

