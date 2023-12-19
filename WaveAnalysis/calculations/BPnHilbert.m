function [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass,varargin)
%BPNHILBERT bandpasses input data according to bandpass and calculates the
%hilbert transform of the BPed data. 
%If there is only one output arg BPnHilbert will only bandpass
%   Usage: 
%   For just BP: FD = BPnHilbert(data,bandpass);
%   For BP+hilbert:  [FD,HT,HTabs,HTangle] = BPnHilbert(data,bandpass);
%   where FD (chXtrialXsample) is the filtered data, HT is the complex hilbert analytic,
%   HTabs is the complex magnitude of the analytic (the envelope) and HT
%   angle is the hilbert phase.
%   For Low/High pass use bandpass=[0 LP] or bandpass=[HP 0] respectively
%   Varargin:
%       -   usePassStop - Use pass and stop frequencies for BP
%       frequencies. When False (default), filter is designed using 
%       highPassCutoff,lowPassCutoff. When true it uses lowPassStopCutoff, 
%       lowPassPassCutoff, highPassStopCutoff and highPassPassCutoff/
%       -   cutwidths: when usePassStop is true, 
%       F.highPassStopCutoff=bandpass(1)-cutwidths(1) and
%       lowPassStopCutoff=bandpass(2)+cutwidths(2). Default is [2 2].
%       -   SamplingFrequency - self explanatory
%       - smoothFD - if true, spatially smoothes FD before calculating
%       Hilbert transform (only for non-average for now). Default is false.
%       requires En.
%       - En - needed for smoothing - electrode map
%       - GaussSigma - gaussian width. default 3

%TODO: QA high pass

usePassStop=0;
cutwidths=[2 2]; %
SamplingFrequency=20000;
smoothFD=false;
GaussSigma=3;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

F=filterData(SamplingFrequency);
F.padding=true;
if all(bandpass) %no zeros
    if usePassStop
        F.lowPassStopCutoff=bandpass(2)+cutwidths(2); %low-pass cutoff frequency
        F.lowPassPassCutoff=bandpass(2);
        F.highPassStopCutoff=bandpass(1)-cutwidths(1); %high-pass cutoff frequency
        F.highPassPassCutoff=bandpass(1);
    else
        F.highPassCutoff=bandpass(1);
        F.lowPassCutoff=bandpass(2);
    end
    F=F.designBandPass;
elseif bandpass(1)==0 %low pass filter
    F.lowPassCutoff=bandpass(2);
    F=F.designLowPass;
else
    F.highPassCutoff=bandpass(1);
    F=F.designHighPass;        
end

FD=F.getFilteredData(data);

if smoothFD
    hsize=size(En);
    h=fspecial('gaussian',hsize,GaussSigma);
    for i=1:size(FD,2)
        FDmovie=convertChannelsToMovie(squeeze(FD(:,i,:)),En,'BGVal',0,'flipEn',0);
        FDmovieFiltered=imfilter(FDmovie,h);
        FD(:,i,:)=convertMovieToChannels(FDmovieFiltered,En);
    end
end

if nargout>1
    HT=zeros(size(FD));
    for i=1:size(FD,2)
        HT(:,i,:)=hilbert(squeeze(FD(:,i,:))').';
    end

    HTabs=abs(HT);
    HTangle=angle(HT);
end
end

