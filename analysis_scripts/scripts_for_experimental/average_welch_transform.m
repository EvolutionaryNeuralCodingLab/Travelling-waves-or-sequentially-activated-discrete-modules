clear all
close all

%% load final results for figure

load([get_wave_analysis_code_base_path() 'precalculated_mats/avg_welch_transform.mat'],'freqs','meanPxxs','resolutonReqd','nTrials')

%% Plot avg welch transform

%plot average <10Hz freqs
plot(freqs(1:round(10/resolutonReqd)),10*log10(meanPxxs(1:round(10/resolutonReqd))))
title(['Trials 1:' num2str(nTrials) ' - AVG Welch Transform'])
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')


%% Calc avg welch transform

data_paths

recObj=binaryRecording(path_to_U4_recording);
load('layout_100_12x12.mat','En')
triggers=recObj.getTrigger;
trigs=tTrig{5};

window_ms=2000;
widenBy=2000; %ms
windenBySamples=widenBy*recObj.samplingFrequency/1000;
band=[0 2];

trialsInBatch=100;
nBatches=10;
nTrials=nBatches*trialsInBatch;

pxxsSum=zeros(100001,recObj.totalChannels,trialsInBatch);

specDelay=250; %ms. only calculate spectrogram starting specDelay after trial start
segmentSize=round((window_ms*recObj.samplingFrequency/1000)/2);
segmentSize_ms=segmentSize/recObj.samplingFrequency*1000;
overlap=round(0.5*segmentSize);
fs=recObj.samplingFrequency;
resolutonReqd = 0.1; %Hz
NFFT = fs / resolutonReqd;


for batch=1:nBatches
    batch
    trials=(1:trialsInBatch)+(batch-1)*trialsInBatch;
    startTimes=trigs(trials);
    [data,time]=recObj.getData([],startTimes+specDelay,window_ms);
    freqs=zeros(100001,1);
    pxxs=zeros(100001,recObj.totalChannels,trialsInBatch);
    for i=1:trialsInBatch
        [pxxs(:,:,i),freqs(:,1)] = pwelch(squeeze(data(:,i,:))',segmentSize,overlap,NFFT,fs);
    end
    pxxsSum=pxxsSum+pxxs;
end
meanPxxs=squeeze(mean(squeeze(mean(pxxsSum,2)),2))/nBatches;

save([get_wave_analysis_code_base_path() 'precalculated_mats/avg_welch_transform.mat'],'freqs','meanPxxs','resolutonReqd','nTrials')
