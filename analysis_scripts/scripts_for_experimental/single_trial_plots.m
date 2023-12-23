clear all
close all

%% Plot signals - raw, HP, LP, Hilbert & ALSA

linesWidths=1;
scaleBarWidth=0.5;

load([get_wave_analysis_code_base_path() 'precalculated_mats/single_trial_plots_data_signals.mat'],...
    'time','single_ch_data','single_ch_FD','single_ch_HT_angle','single_ch_HP','single_ch_bin_spikes', ...
    'samplingFrequency','single_ch_ALSA'...
    )


%%%
f=figure;
h(1)=plot(time,single_ch_data,'Color',0.7*[1 1 1],'LineWidth',linesWidths);
hold on
h(2)=plot(time,single_ch_FD,'b','LineWidth',linesWidths);
plot([50,50],[100,150],'b','LineWidth',scaleBarWidth) %scale bar
add2phase=-abs(min(single_ch_data))-max(single_ch_HT_angle*10)-5;
h(3)=plot(time,single_ch_HT_angle*10+add2phase,'Color',[0 0.9 0],'LineWidth',linesWidths);
add2scalebar=add2phase-abs(min(single_ch_HT_angle)*10);
plot([50,50],[add2scalebar+20*pi,add2scalebar+20*pi+2*pi*10],'Color',[0 0.9 0],'LineWidth',scaleBarWidth) %scale bar
add2HP=+max(single_ch_data)+abs(min(single_ch_HP)+10);
h(4)=plot(time,squeeze(single_ch_HP)+add2HP,'k','LineWidth',linesWidths);
singleChannelSpikes=find(single_ch_bin_spikes)/samplingFrequency*1000;
add2spikes=add2HP+max(single_ch_HP)+10;

for i=1:length(singleChannelSpikes)
   plot([singleChannelSpikes(i) singleChannelSpikes(i)],[add2spikes add2spikes+15],'k','LineWidth',0.7)
end
add2alsa = add2spikes + 50;
h(5) = plot(time, single_ch_ALSA+add2alsa,'r','LineWidth',linesWidths);
xlabel('Time [ms]')
leg=legend(h,{'Unfiltered Data','LP Filtered Data','Hilbert Phase','HP Filtered Data','ALSA'},'Location','north','Orientation','horizontal','NumColumns',5);
legend('boxoff')
ylimit=ylim;
ylim([ylimit(1) add2alsa+155])
yticks([])
%manually set legend size cuz its huge
ax=gca;
ax.Legend.FontSize=7.5;
ax.Legend.ItemTokenSize=[5 48];


%% Plot Latency maps - LPF phase crossings and ALSA

load([get_wave_analysis_code_base_path() 'precalculated_mats/single_trial_plots_data_latency_maps.mat'],...
    'chs','normedTimes_LFP','En','waveCenterPath','relevantChannels','normedTimes_ALSA')

%Plot LFP PLM
f=figure;
[hCbar,h]=IntensityPhysicalSpacePlot(chs,normedTimes_LFP,En,'plotElectrodeNumbers',0,'plotSizeBar',0,'plotGridLines',0,'markerSize',10);
set(h,'Units','centimeters','Position',[1.5,0.5,1.5,1.5])
set(hCbar,'Units','centimeters','Position',[1.2,0.5,0.15,1.5])

ylabel(hCbar,{'Latency','[ms]'},'Units','centimeters','Position',[-0.25,0.71,0],'FontSize',8);
set(f,'Units','centimeters','Position',[18 23.3 3.3 2.4]);
hold on
scatter(waveCenterPath(:,1),waveCenterPath(:,2),10,'k')

%Plot ALSA PLM
f=figure;
[hCbar,h]=IntensityPhysicalSpacePlot(relevantChannels,normedTimes_ALSA,En,'plotElectrodeNumbers',0,'plotSizeBar',0,'plotGridLines',0,'markerSize',10);
set(h,'Units','centimeters','Position',[1.5,0.5,1.5,1.5])
set(hCbar,'Units','centimeters','Position',[1.2,0.5,0.15,1.5])

ylabel(hCbar,{'Latency','[ms]'},'Units','centimeters','Position',[-0.25,0.71,0],'FontSize',8);
set(f,'Units','centimeters','Position',[18 23.3 3.3 2.4]);
hold on
scatter(waveCenterPath(:,1),waveCenterPath(:,2),10,'k')


%% U4 trial 17 calculations - load for all calcualtions

data_paths

% Choose trial and load

trial=17;

recObj=binaryRecording(path_to_U4_recording);
load('layout_100_12x12.mat','En')
nCh=max(En(:));

tTrig=recObj.getTrigger;
trigs=tTrig{5};

startTimes=trigs(trial);

%% Calculate signals - raw, HP, LP, Hilbert & ALSA

window_ms=3000;
widenBy=2000; %ms
widenBySamples=widenBy*recObj.samplingFrequency/1000;

binSpikes = getSpikeBinMatByChannel(path_to_U4_tIc,startTimes(1),startTimes(1)+window_ms,recObj.samplingFrequency);

band=[0 2];
HPband=[200 0];

[data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
[~,~,FD_HP,HT_HP,HTabs_HP,HTangle_HP] = getCroppedFD(recObj,startTimes,window_ms,widenBy,HPband);
spikeRateSmoothing=2000;
nNearestChannels=4;
ALSA = getALSAFromTIC(path_to_U4_tIc,startTimes,window_ms,En,recObj.samplingFrequency,'spikeRateSmoothing',spikeRateSmoothing,'nNearestChannels',nNearestChannels);

%   plot raw, FD and hilbert   %
singleChannel=100; 

single_ch_data = squeeze(data(singleChannel,1,:));
single_ch_FD = squeeze(FD(singleChannel,1,:));
single_ch_HT_angle = squeeze(HTangle(singleChannel,1,:));
single_ch_HP = squeeze(FD_HP(singleChannel,1,:));
single_ch_bin_spikes = binSpikes(singleChannel,:);
samplingFrequency = recObj.samplingFrequency;
single_ch_ALSA = ALSA(singleChannel,:);


save([get_wave_analysis_code_base_path() 'precalculated_mats/single_trial_plots_data_signals.mat'],...
    'singleChannel','time','single_ch_data','single_ch_FD','single_ch_HT_angle','single_ch_HP','single_ch_bin_spikes', ...
    'samplingFrequency','single_ch_ALSA'...
    )

%% Latency maps - LPF phase crossings and ALSA

% defs and loads
window_ms=1500;
widenBy=2000; %ms
ms2samples=recObj.samplingFrequency/1000;
widenBySamples=widenBy*ms2samples;

band=[0 2];

excitationAngle=90;
excitationPhase=excitationAngle*pi/180;
crossingType=4;
crossingsName={'Maxima','Minima','Halfway Up','Excitaion'};

maxTempDist=2000;
minChannelInWave=80;
minHilbertAmp=10;

nCh=max(En(:));
chPos=calcChannelsPosition(En);

% calcs
[data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band);
[crossings,hilbertAmps] = getHilbertCrossings(squeeze(HTabs(:,1,:)),squeeze(HTangle(:,1,:)),'excitationPhase',excitationPhase);
binSpikes = getSpikeBinMatByChannel(path_to_U4_tIc,startTimes(1),startTimes(1)+window_ms,recObj.samplingFrequency);
[clusterLimits,channels,times] = getTrialClusters(crossings{crossingType},En,maxTempDist,minChannelInWave,binSpikes,'plotTrialsClusters',0,'hilbertAmps',hilbertAmps{crossingType},'minHilbertAmp',minHilbertAmp);
startEndWave=clusterLimits(1,:);
waveWidth=startEndWave(2)-startEndWave(1);
ALSALimits=[max(1,startEndWave(1)-waveWidth/2),min(round(window_ms*recObj.samplingFrequency/1000),startEndWave(2)+waveWidth/2)];
[ALSA_Locs,channelsWithALSA] = getALSAPeaksFromTIC(path_to_U4_tIc,startTimes(1),window_ms,En,recObj.samplingFrequency,'startEndWave',ALSALimits,'onsetType','firstMax'); 
[relevantChannels,relevantCrossingTimes,relevantALSATimes]=getRelevantWaveTimes(channels{1},times{1},channelsWithALSA,ALSA_Locs,nCh);
waveCenterPath = drawWavePath(crossings{crossingType},hilbertAmps{crossingType},startEndWave,En,'normCoordinates',0,'flipEn',0);

chs = channels{1};
normedTimes_LFP=(times{1}-min(times{1}))/recObj.samplingFrequency*1000; %when possible, should be relevantCrossingTimes from getRelevantWaveTimes. Also channels{1} should be relevantChannels
normedTimes_ALSA=(relevantALSATimes-min(relevantALSATimes))/recObj.samplingFrequency*1000;

save([get_wave_analysis_code_base_path() 'precalculated_mats/single_trial_plots_data_latency_maps.mat'],...
    'chs','normedTimes_LFP','En','waveCenterPath','relevantChannels','normedTimes_ALSA')