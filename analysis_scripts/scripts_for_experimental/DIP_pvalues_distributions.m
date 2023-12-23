%% DIP-test p-values - LPF vs ALSA

clear all
close all

%% Load final results for figure

load([get_wave_analysis_code_base_path() 'precalculated_mats/U4Statistics90phase.mat'],'dip_pvalues','trialsParticipated')

%% % Plot dip p-values

[f,PoD]=plotPvalues(dip_pvalues(1,trialsParticipated),dip_pvalues(2,trialsParticipated),'plotLines',0);
set(gcf,'Units','centimeters','Position',[14.0000   12.9381    9.3    5.45])

%% Calculation of DIP-test p-values
data_paths
recObj=binaryRecording(path_to_U4_recording);

triggers=recObj.getTrigger;
trigs=triggers{5};

load('layout_100_12x12.mat','En')
nCh=max(En(:));
chPos=calcChannelsPosition(En);

window_ms=1500;
widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

band=[0 2];

nTrials=length(trigs);
nTrialsInBatch=100;
nBatches=ceil(nTrials/nTrialsInBatch);
allTrials=1:nTrials;

trialsParticipated=false(1,nTrials);
dip_pvalues=zeros(2,nTrials);

crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitaion'};
   
excitationPhase=90*pi/180; 


maxTempDist=2000;
minChannelInWave=80;
minHilbertAmp=10;

for batch=1:nBatches
    disp(num2str(batch))
    
    trials=allTrials(((batch-1)*nTrialsInBatch+1):min(batch*nTrialsInBatch,nTrials));
     startTimes=trigs(trials);
    
    [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
    
    for i=1:length(trials)
        trial=trials(i);
        [crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,i,:)),squeeze(HTangle(:,i,:)),'excitationPhase',excitationPhase);
        binSpikes = getSpikeBinMatByChannel(path_to_U4_tIc,startTimes(i),startTimes(i)+window_ms,recObj.samplingFrequency);
        if sum(binSpikes(:))<120 %if there is less than 1 spike per ch on average, skip this trial
            continue
        end
        [clusterLimits,channels,times,spikesPerCluster,allSeedSamples,allSeedChannels] = getTrialClusters(crossings{crossingType},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'minHilbertAmp',minHilbertAmp);
        if isempty(clusterLimits)
            continue
        end
        startEndWave=clusterLimits(1,:);
        waveWidth=startEndWave(2)-startEndWave(1);
        if startEndWave(1)>round(window_ms*recObj.samplingFrequency/1000)/2 %if the pattern starts after half of the trial, this is probably false pattern
            continue
        end
       
        ALSALimits=[max(1,startEndWave(1)-waveWidth/2),min(round(window_ms*recObj.samplingFrequency/1000),startEndWave(2)+waveWidth/2)];
        [ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(path_to_U4_tIc,startTimes(i),window_ms,En,recObj.samplingFrequency,'startEndWave',ALSALimits,'onsetType','firstMax');
        [relevantChannels,relevantCrossingTimes,relevantALSATimes]=getRelevantWaveTimes(channels{1},times{1},channelsWithALSA,ALSA_Locs,nCh);

        if length(relevantChannels)<4
            continue
        end
        waveCenterPath = drawWavePath(crossings{crossingType},hilbertAmps{crossingType},startEndWave,En,'normCoordinates',1,'flipEn',0);
        nSamples=size(waveCenterPath,1);
        chPosMax=max(chPos(:));

        ALSACoordinates=[chPos(relevantChannels,2)/chPosMax chPos(relevantChannels,1)/chPosMax (relevantALSATimes-min(relevantALSATimes))'/(max(relevantALSATimes)-min(relevantALSATimes))];
        crossingsTimes=(relevantCrossingTimes-min(relevantCrossingTimes))'/(max(relevantCrossingTimes)-min(relevantCrossingTimes));
        crossingsCoordinates=[chPos(relevantChannels,2)/chPosMax chPos(relevantChannels,1)/chPosMax crossingsTimes];


        [~, p_LFP_temporal] = hartigansdipsigniftest(sort(relevantCrossingTimes), 500);
        [~, p_ALSA_temporal] = hartigansdipsigniftest(sort(relevantALSATimes), 500);

        dip_pvalues(1:2,(batch-1)*nTrialsInBatch+i)=[p_LFP_temporal;p_ALSA_temporal];
        trialsParticipated((batch-1)*nTrialsInBatch+i)=true;
    end
   save([get_wave_analysis_code_base_path() 'precalculated_mats/tempU4Statistics90phase.mat'],'dip_pvalues','trialsParticipated')    
end

save([get_wave_analysis_code_base_path() 'precalculated_mats/U4Statistics90phase.mat'],'dip_pvalues','trialsParticipated')