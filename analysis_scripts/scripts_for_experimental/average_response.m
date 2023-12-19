% addpath(genpath('/media/E/Yuval/YuvalNET2'));

clear all
close all

data_paths

%% Calc (or simply load in next section) avg response
%Calc
recObj=binaryRecording(path_to_U4_recording);
load('layout_100_12x12.mat','En')
load(path_to_U4_digital_triggers,'tTrig')
trigs=tTrig{5};

window_ms=1500;
widenBy=2000; %ms
windenBySamples=widenBy*recObj.samplingFrequency/1000;
band=[0 2];

trials=1:1000;
trialsInBatch = 100;
nTrials=length(trials);
startTimes=trigs(trials);


[~,time,FDmean,HTmean,HTabsmean,HTanglemean] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band,'returnAVG',1,'trialsInBatch',trialsInBatch);
save([get_wave_analysis_code_base_path() 'precalculated_mats/avg_response.mat'],'time','FDmean','HTmean','HTabsmean','HTanglemean')

%% Load results from last section
load([get_wave_analysis_code_base_path() 'precalculated_mats/avg_response.mat'],'time','FDmean','HTmean','HTabsmean','HTanglemean')

%% Find peaks and plot
load('layout_100_12x12.mat','En')
[crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabsmean(:,1,:)),squeeze(HTanglemean(:,1,:)));
crossingType=1;
[clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,2000,80,[],'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType});
startEndWave=clusterLimits(1,:);
waveWidth=startEndWave(2)-startEndWave(1);

[relevantChannels,relevantCrossingTimes]=getClusterFirstCrossings(channels{1},times{1});

waveCenterPath = drawWavePath(crossings{crossingType},hilbertAmps{crossingType},startEndWave,En,'normCoordinates',0);
[roundWaveCenterPath,centerMoveTimes,pathChannels,waveChannelsPos] = linearizeWaveCenterPath(waveCenterPath,En);

%find actual maxima (not the phase crossing)
nRelevant=length(relevantCrossingTimes);
relevantMaximaTimes=zeros(1,nRelevant);
for i=1:nRelevant
    [pks,locs] = findpeaks(squeeze(FDmean(relevantChannels(i),1,:)));
    [~,maxInd]=max(pks);
    relevantMaximaTimes(i)=locs(maxInd);
end

samplingFrequency=20e3;
time=(1:size(FDmean,3))/samplingFrequency*1000; %sampling rate of 20kHz
plotWaveFD(time,FDmean(:,1,:),relevantChannels,relevantMaximaTimes,pathChannels,samplingFrequency,'scatterSize',45);

% plot filmstrip with scatter on first frame
flimstripStartEnd=[200 750]*samplingFrequency/1000;
plotFilmstrip(FDmean(:,1,:),flimstripStartEnd,En,waveChannelsPos,samplingFrequency,'colorbarFontSize',8); 

