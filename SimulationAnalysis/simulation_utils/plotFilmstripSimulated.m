function [f] = plotFilmstripSimulated(data,flimstripStartEnd,En,crossings,amps,varargin)
%PLOTFILMSTRIP plots a filmstrip of the data in data, at fixed intervals
%inside filmestripStartEnd
%   INPUT:
%   - data (EnHeightXEnWidthXnSamples) - the data to draw (this is
%   different format than plotFilmstrip
%   - filmstripStartEnd (1X2) - start and end samples for the
%   filmstrip. 8 frames equally seperated are taken from this.
%   - En - electrode layout by which the frames are set up from data
%   - waveChannelsPos (nElectrodesX2) - x,y positions of electrodes to
%   plot. flipped ud.
%   - varargins
%       nRows - number of rows in subplot. Default is 1
%       nTimes - number of frames to plot. Default is 8

nRows=1;
nTimes=8;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if ndims(data)~=3 || size(data,1)~=size(En,1)|| size(data,2)~=size(En,2)
   warning('Problem with data array. Should be in movie format!') 
end

f=figure;

representativeTimes=round(linspace(flimstripStartEnd(1),flimstripStartEnd(2),nTimes)); %8 snippets
nRep=length(representativeTimes);
minMovieVal=min(min(min(data(:,:,representativeTimes))));
maxMovieVal=max(max(max(data(:,:,representativeTimes))));
% ax1=subplot(nRows,ceil(nRep/nRows),1);
% imagesc(ax1,data(:,:,representativeTimes(1)),[minMovieVal maxMovieVal])
% 
% 
% P = get(ax1,'Position');
% XLIM = get(ax1,'XLim');
% YLIM = get(ax1,'YLim');
% PA = get(ax1,'PlotBoxAspectRatio');
% set(ax2,'Position',P,'XLim',XLIM,'YLim',YLIM,'PlotBoxAspectRatio',PA)
% 
% 
% %%Link them together
% linkaxes([ax1,ax2])
% %%Hide the top axes
% ax1.Visible = 'off';
% ax1.XTick = [];
% ax1.YTick = [];
% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];
% %%Give each one its own colormap
% colormap(ax1,'parula')
% colormap(ax2,'jet')
% %%Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
% cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
% title(['T=' num2str(round(representativeTimes(1)/samplingFrequency*1000)) 'ms'])

for i=1:nRep
    subplot(1,ceil(nRep/nRows),i)
    imagesc(data(:,:,representativeTimes(i)),[minMovieVal maxMovieVal])
%     title(['T=' num2str(round(representativeTimes(i)/samplingFrequency*1000)) 'ms'])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
end

% hp4 = get( subplot(2,ceil(nRep/2),i),'Position');
% colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*2.1])

%     pos=get(gcf,'Position')
set(gcf,'Position',[300   862   981   105]);


end

