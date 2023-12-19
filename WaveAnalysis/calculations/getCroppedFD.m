function [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,startTimes,window_ms,widenBy,band,varargin)
% Example usage: [data,time,FD,HT,HTabs,HTangle] = getCroppedFD(recObj,trig(i),1500,2000,[0 2])
%GETCROPPEDFD BP filters and performs hilbert transform on data that starts at
%startTimes (ms) for the duration of window_ms. It does so by applying the
%filter and hilbert on a larger temporal window, widenBy (ms) before and
%after the start and end time points of the data, and then cropping.
%This mean that to filter the data from 1s to 1.5s it first filters
%1-widenBy to 1.5+widenBy and then croppes the padding data.
%   INPUT:
%   - recObj - recording object
%   - startTimes - vector of times (ms) from the start of the recording
%   - window_ms - temporal window we want to filter/hilbert
%   - widenBy - (ms) how much paddin to do
%   - band - 1X2 array - BP range (Hz) for Lowpass use e.g. [0 2].
%   - Possible Varargins (given as 'Key',Value pairs)
%       - returnAVG (logical) - If false (default), performs the hilbert 
%       transform and crossing analysis on each FD trial. If true, it first
%       averages the FD and performs analysis on average. 
%       Notice - this will return the average FD, but the full data array
%       (i.e. data is nChannelsXnTrialsXnSamples while FD is nChannelsX1XnSamples
%       - trialsInBatch - if returnAVG is set to 1, this variable can be
%       used to calculate the average in batches. Default is all trials,
%       i.e. length(startTimes). This must be a divisor of numer of trials.
%       - smoothFD - if true, spatially smoothes FD before calculating
%       Hilbert transform (only for non-average for now). Default is false.
%   OUTPUT:
%       if returnAVG=1, data will only contain the data from last batch
%
% NOTICE: Currently this function gets data for all channels. If at some
% point the option to choose only selected channels will be implemented,
% make sure to notice that this may cause problem when applying hilbert
% transform on squeezed 1X1XnSamples arrays (i.e. of single channels), 
% because it will create an nSamplesX1 array - so make sure all transposes 
% are needed). Also, if at some point function will be change to accommade
% trialsInBatch the is not a divisor of number of trials, make sure to
% no errors in the averaging process (e.g. deviding each batch by the
% number of trials in it, includin last batch, and also making sure that
% empty trials in last batch are not averaged in "mean(FD,2)"

nTrials=length(startTimes);
returnAVG=0;
trialsInBatch=nTrials;
smoothFD=false;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if mod(nTrials,trialsInBatch)>0
   error('trialsInBatch is not a divisor of the number of trials')
end
nBatches=ceil(nTrials/trialsInBatch);

window_sample=window_ms*recObj.samplingFrequency/1000;
windenBySamples=widenBy*recObj.samplingFrequency/1000;

if ~returnAVG
    [data,time]=recObj.getData([],startTimes-widenBy,window_ms+2*widenBy);
    if smoothFD
        if exist('GaussSigma','var')
            [FD,HT,HTabs,HTangle] = BPnHilbert(data,band,'smoothFD',smoothFD,'En',En,'GaussSigma',GaussSigma);
        else
            [FD,HT,HTabs,HTangle] = BPnHilbert(data,band,'smoothFD',smoothFD,'En',En);
        end
    else
         [FD,HT,HTabs,HTangle] = BPnHilbert(data,band);
    end
else
    FDsum=zeros(recObj.totalChannels,trialsInBatch,window_sample+2*windenBySamples);
    fprintf('\nBatch Num (/%d) : ',nBatches);
    nDigits=0;
    for batch=1:nBatches
        fprintf([repmat('\b',1,nDigits) '%d'],batch);nDigits=length(num2str(batch));
        trials=(1:trialsInBatch)+trialsInBatch*(batch-1);
        batchStartTimes=startTimes(trials);    
        [data,time]=recObj.getData([],batchStartTimes-widenBy,window_ms+2*widenBy);
        FDbatch = BPnHilbert(data,band);
        FDsum=FDsum+FDbatch;
    end
    fprintf('\n')
    FD = mean(FDsum,2)/nBatches;
    HT(:,1,:)=hilbert(squeeze(FD(:,1,:))').';
    HTabs=abs(HT);
    HTangle=angle(HT);
end
%crop spare data
data=data(:,:,windenBySamples+1:end-windenBySamples);
time=time(1:end-2*windenBySamples);
FD=FD(:,:,windenBySamples+1:end-windenBySamples);
HT=HT(:,:,windenBySamples+1:end-windenBySamples);
HTabs=HTabs(:,:,windenBySamples+1:end-windenBySamples);
HTangle=HTangle(:,:,windenBySamples+1:end-windenBySamples);

end

