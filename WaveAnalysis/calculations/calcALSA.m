function ALSA=calcALSA(spikingRate,varargin)
%CALCALSA calculates the average local spiking activity (ALSA) for every
%channel. For every channel ch and for every sample s ALSA is defined by
%the weighted average of nNearestChannels. If En is given, nearest
%neighbors are locally adjacent channels and nNearestChannels must be
%either 4 (4 nearest neighbors,default) or 8 (4 nearest and 4 next nearest
%neighbors). If En is not given, the weighted average is calculated using
%nChs2average (default: 2) above and below channel ch (with cyclic channel counting) 

%   INPUT:
%       - spikingRate (nChXnSamples) - spiking rate to average spatially.
%       - Varargs (given as 'Name','Value' pairs):
%           - nChs2average - if En is not given as varargin, the weighted 
%           average is calculated using nChs2average above and below 
%           channel ch (with cyclic channel counting). Default is 2.
%           - Weights: if En is not given as varargin, this will be the
%           vector of averaging weights. Must be of length 
%           nChs2average*2+1. Default is gausswin(nChs2average*2+1)/sum(gausswin(nChs2average*2+1))
%           - En - Channel Layout. If given averaging will be over nearest 
%           neighbors.
%           - nNearestChannels - number of channels to average over. Can be
%           either 4 (default) or 8 (4nn+4nnn). weights are 1 for the
%           center channel, 1/2 for nn and 1/4 for nnn (normalized)
%   OUTPUT:
%       - ALSA (nChXnSamples)

nChs2average=2;
Weights=gausswin(nChs2average*2+1)/sum(gausswin(nChs2average*2+1));
nNearestChannels=4; 

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end
    
[nCh,nSamples]=size(spikingRate);

ALSA=zeros(nCh,nSamples);

if exist('En','var') 
    if nNearestChannels==4
        for i=1:nCh
            [nn,~]=getNeighbors(i,En);
            ALSA(i,:)=(spikingRate(i,:)+sum(spikingRate(nn,:))*0.5)/(1+0.5*length(nn)); %give nn half the weight, and average by total weight (1+0.5*4)
        end
    elseif nNearestChannels==8
         for i=1:nCh
            [nn,nnn]=getNeighbors(i,En);
            ALSA(i,:)=(spikingRate(i,:)+sum(spikingRate(nn,:))*0.5+sum(spikingRate(nnn,:))*0.25)/(1+0.5*length(nn)+0.25*length(nnn)); %give nn half the weight,nnn quarter weight, and average by total weight (1+0.5*4+0.25*4)
        end
    else
        error('If Electrode En array is given, nNearestChannels must be either 4 or 8')
    end
else %no En given, calculate average according to numerical proximity

    if size(Weights,2)>1
        error('Weights must be a column vector!')
    end

    if length(Weights)~=(2*nChs2average+1)
        error('Weights vector must be 2*nChs2average+1 long!')
    end

    if abs(sum(Weights)-1)>1e-15 %up to rounding error
        warning('Weights vector is not normalized. I''ll do it myself')
        Weights=Weights/sum(Weights);
    end

    for i=1:nCh
        ch2average=1+mod(((i-nChs2average):(i+nChs2average))-1,nCh);
        ALSA(i,:)=sum(spikingRate(ch2average,:).*Weights);
    end
end

end

