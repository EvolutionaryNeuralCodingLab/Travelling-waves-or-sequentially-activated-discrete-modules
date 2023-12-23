%% correlations statistics and shuffling

clear all
close all

%% Load final results for figure

load([get_wave_analysis_code_base_path() 'precalculated_mats/correlations.mat'],'rCoefs','rShuffles')

%% Plot PLDC vs shuffled

edges=linspace(-1,1,31);
histogram(rCoefs,edges,'Normalization','Probability');
hold on
histogram(rShuffles(:),edges,'Normalization','Probability');

xlabel('PLDC')
ylabel('Frequency')
legend({'Real data','Shuffled Data'},'Location','northwest')
legend('box','off')
ax=gca;
ax.Legend.FontSize=6;
ax.Legend.ItemTokenSize=[5 48];

ylim([0 0.6])
set(gcf,'Units','centimeters','Position',[17 14 4 3]);

%% Calculation of PLDC vs shuffled

data_paths

window_ms=1500; %ms
widenBy=2000;
nCh=120; %number of channels - in code this will channels arrays will be 1:nCh
band=[0 2];


recObj=binaryRecording(path_to_U4_recording);

triggers=recObj.getTrigger;
trigs=triggers{5};

load('layout_100_12x12.mat','En')
excitationPhase=90*pi/180;
crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitation'};

trialsInBatch=100;
nBatches=10;
nTrials=nBatches*trialsInBatch;

nShuffles=100;
rCoefs=nan(1,nTrials);
distances=cell(1,nTrials);
allTrialsTimes=cell(1,nTrials);
allTrialsChannels=cell(1,nTrials);
shuffledRMean=nan(1,nTrials);
shuffledRstd=nan(1,nTrials);
rShuffles=nan(nShuffles,nTrials);

for batch=1:nBatches
    batch
    trials=(1:trialsInBatch)+(batch-1)*trialsInBatch;
    startTimes=trigs(trials);
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);  
    for i=1:trialsInBatch
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        binSpikes = getSpikeBinMatByChannel(path_to_U4_tIc,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency,nCh);
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,1500,80,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'plotStyles',{'b.','g.'});
        if ~isempty(clusterLimits)
            [rCoefs(trials(i)),distances{trials(i)}] = calcDistanceAndPhaseLatencyCorrelation(channels{1},times{1},En);
            allTrialsTimes{trials(i)}=times{1};
            allTrialsChannels{trials(i)}=channels{1};
            for j=1:nShuffles
                [rShuffles(j,trials(i)),~] = calcDistanceAndPhaseLatencyCorrelation(channels{1},times{1}(randperm(length(times{1}))) ,En);
            end
            shuffledRMean(i)=mean(rShuffles(:,trials(i)));
            shuffledRstd(i)=std(rShuffles(:,trials(i)));
        end
    end
    save([get_wave_analysis_code_base_path() 'precalculated_mats/correlations.mat'],'rCoefs','allTrialsTimes','allTrialsChannels','distances','shuffledRMean','shuffledRstd','rShuffles')

end