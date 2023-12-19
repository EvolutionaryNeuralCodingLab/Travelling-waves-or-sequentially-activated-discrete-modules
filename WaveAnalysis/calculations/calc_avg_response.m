function [time,FDmean,HTmean,HTabsmean,HTanglemean] = calc_avg_response(recObj,digital_trig_path)
%CALC_AVG_RESPONSE Summary of this function goes here
%   Detailed explanation goes here

tTrig = load(digital_trig_path,'tTrig');
trigs=tTrig{5};

window_ms=1500;
widenBy=2000; %ms
windenBySamples=widenBy*recObj.samplingFrequency/1000;
band=[0 2];

% export average 

% trials=1:1000;
% trialsInBatch = 100;
trials=1:10;
trialsInBatch = 10;
nTrials=length(trials);
startTimes=trigs(trials);


[~,time,FDmean,HTmean,HTabsmean,HTanglemean] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band,'returnAVG',1,'trialsInBatch',trialsInBatch);


end

