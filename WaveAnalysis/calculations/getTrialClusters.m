function [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings,En,maxTempDist,minChannelInWave,binSpikes,varargin)
%FINDCONTINUOUSCLUSTERS finds the start and end times of spatiotemporal 
% continuous crossing clusters, representing waves. 
%       Example usage:
%   [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,1500,80,binSpikes,'plotTrialsClusters',1,'hilbertAmps',hilbertAmps{crossingType},'plotStyles',{'b.','g.'});

%       Input:
%       -   crossings: nChXnCrossings crossings times (one of the 4 
%       "crossings" matrices given by getHilbertCrossings)
%       -   En: channel layout 
%       -   maxTempDist: the maximal time window allowed for two neighbors 
%       to have crossings
%       -   minHilbertAmp - Only look at crossings above this threshold.
%       Default it 0 (no threshold)
%       -   binSpikes - nChXnSamples logical matrix with ones marking spike 
%           times. This is the output of getSpikeBinMatByChannel function.
%           If empty, spikesPerCluster will have be all zeros
%       -   Varargins (given as 'Key',Value pairs):
%           -   minChannelInWave: minimal channels to be included in a cluster
%           -   minSpikesPerCluster - If given, function only returns 
%           cluster with more spikes within cluster limits than this.
%           -   plotTrialsClusters (logical) - plot the results of this
%           trial crossing clustering. Fro this hilbertAmps also should be
%           sent as varargin. deafult is false.
%           -   hilbertAmps - nChXnCrossings hilbert amplitudes at crossings the
%           times of the crossings (one of the 4 hilbertAmps matrices given by 
%           "getHilbertCrossings")
%           -   plotStyles - cell array of styles to use for each cluster.
%           Only relevant if plotTrialsClusters=1. If so, cluster i will be
%           plotted using the style plotStyles{((i-1) mod numel(plotStyles))+1}.
%           So if you want triplets of clusters to be bluee,red,green dots
%           repeatedly then use {'.b','.r','.g'}. Default is {'.b'}. 
%           -   plotSpikes - logical - Default is true.
%           -   timeInms - logical. Display crossings times in ms, using
%           sample2ms. Default is 0.
%           -   sample2ms - 1000/samplingFrequency
%           -   plotColorbar - logical
%           -   sz - circles size. default is 25.
%           -   markersize - detected clusters' markersize
%           -   order
%
%       Output:
%       -   clusterLimits (nGoodClustersX2) - position (samples) of the
%       first and last crossings in the cluster.
%       -   channels,times - 1XnGoodClusters cell array with arrays of the
%       channels participated in each cluster and the times of the
%       crossings (channel j of cluster i was channels{i}(j) and crossed 
%       at time times{i}(j) (all in samples relative to the beginning of 
%       the trial).
%       -   spikesPerCluster (nGoodClustersX1) - vector with number of
%       spikes within each cluster
% TODO: Consider moving plotting option into plotSingleHilbertCrossing

plotTrialsClusters=false;
minHilbertAmp=0;
plotStyles={'b.'};
plotSpikes=1;
timeInms=0;
sample2ms=1; %assuming no conversion
plotColorbar=1;
sz=25;
markersize=6;
order=1:max(En(:));

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

[~,chOrder]=sort(order);

if isempty(binSpikes)
    nSamples=max(crossings(:));
    noSpikes=1;
else
    nSamples=size(binSpikes,2);
    noSpikes=0;
end

%Get rid of crossings with amps lower than minHilbertAmp
if minHilbertAmp>0
    if ~exist('hilbertAmps','var')
        error('In order to use minHilbertAmp threshold, hilbertAmps must be given as varargin. See Function Description')
    else
        lowAmpCrossingsInd=find(hilbertAmps<minHilbertAmp);
        hilbertAmps(lowAmpCrossingsInd)=0;
        crossings(lowAmpCrossingsInd)=0;
    end
end


%convert crossings to binary
sampleCrossings=getCrossingsBySamples(crossings,[],'nSamples',nSamples);

if plotTrialsClusters
    if ~exist('hilbertAmps','var')
        error('hilbertAmps must be given as varagin in order to plot crossings');
    end
    if plotSpikes
        plotSingleHilbertCrossing(crossings,hilbertAmps,0,'',1,'Spikes',binSpikes,'plotLegend',false,'timeInms',timeInms,'sample2ms',sample2ms,'plotColorbar',plotColorbar,'sz',sz,'order',order);
    else
        plotSingleHilbertCrossing(crossings,hilbertAmps,0,'',1,'plotLegend',false,'timeInms',timeInms,'sample2ms',sample2ms,'plotColorbar',plotColorbar,'sz',sz,'order',order);
    end
end


nGoodClusters=0;
spikesPerCluster=[];
spikesPerSample=sum(binSpikes);
clusterLimits=[];
channels ={};
times={};
allSeedChannels=[];
allSeedSamples=[];
[seedCh,seedSample]=find(sampleCrossings,1); %find the first non zero crossings. find searches rows first so all columns left to what it finds will be zeros
while ~isempty(seedCh)
    allSeedChannels=[allSeedChannels seedCh]; %DEBUG THIS: last bug should be deleted if cluster is not included (e.g. if nChInWave>=minChannelInWave)
    allSeedSamples=[allSeedSamples seedSample];
    %get the wave around the seed
    [nChInWave,clusterChannels,clusterTimes] = countContinousCrossings(crossings,En,maxTempDist,seedCh,seedSample);
    
    
%     %%%QAAAA%%%
%     currentClusterLimits=[min(clusterTimes),max(clusterTimes)];
%     h(1)=plot(currentClusterLimits,[0 0],'LineWidth',2,'Color','k');
%     h(2)=plot(clusterTimes,clusterChannels,'.');
%     nChInWave
%     delete(h(1))
%     delete(h(2))
%     %%%%%%%%%%%
    if nChInWave>=minChannelInWave
        currentClusterLimits=[min(clusterTimes),max(clusterTimes)];
        if ~noSpikes
            currentClusterSpikes=sum(spikesPerSample(currentClusterLimits(1):currentClusterLimits(2)));
        end
        if ~exist('minSpikesPerCluster','var') || (exist('minSpikesPerCluster','var') && ~noSpikes && currentClusterSpikes>=minSpikesPerCluster)
            nGoodClusters=nGoodClusters+1;
            clusterLimits(nGoodClusters,1:2)=currentClusterLimits;
            channels{nGoodClusters}=clusterChannels;
            times{nGoodClusters}=clusterTimes;
            if noSpikes
                spikesPerCluster(nGoodClusters)=0;
            else
                spikesPerCluster(nGoodClusters)=currentClusterSpikes;
            end
        end
    end
    %erase current cluster
    sampleCrossings(sub2ind(size(sampleCrossings),clusterChannels,clusterTimes))=0;
    
    [seedCh,seedSample]=find(sampleCrossings,1);
end

nStyles=numel(plotStyles);
if plotTrialsClusters && nGoodClusters>0
    for i=1:size(clusterLimits,1)
        plot(clusterLimits(i,:)*sample2ms,[0 0],'LineWidth',2,'Color','k')
        plot(times{i}*sample2ms,chOrder(channels{i}),plotStyles{mod(i-1,nStyles)+1},'markerSize',markersize)
    end
end
end

