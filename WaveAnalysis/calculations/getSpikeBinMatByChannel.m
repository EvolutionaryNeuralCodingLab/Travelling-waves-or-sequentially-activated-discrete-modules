function binSpikes = getSpikeBinMatByChannel(ticPath,start_times_ms,end_times_ms,sampling_rate,nCh)
%GETSPIKEBINMAT returns a nChXnSamples logical matrix with ones marking
%spike times
%   This function gets the spike times in t,ic format. It goes through all
%   start_times and finds all spikes in time window from start_times(i) to 
%   end_times(i) and marks ones in all the times(/samples) where each 
%   channel fired. 
%   If nCh is not given, it is taken to be the last channel in ic (i.e. ic(1,end))
%   All trial spikes are concatenated to get a sequence in accordance to
%   the matrices of getDataSequence.
%   sampling_rate is samples/s.

load(ticPath,'t','ic');

if ~exist('nCh','var')
    nCh=ic(1,end);
end
nTrials=numel(start_times_ms);

ms2sample=sampling_rate/1000;
window_ms=end_times_ms(1)-start_times_ms(1);
window_samples=round(window_ms*ms2sample);
nSamples=window_samples*nTrials;

binSpikes=false(nCh,nSamples);


for i=unique(ic(1,:))
    iChannelColumns=find(ic(1,:)==i);
%     allChannelSpikes=[];
    for j=iChannelColumns 
        neuronSpikes=t(ic(3,j):ic(4,j));
        for k=1:nTrials
            neuronSpikesInWindow=neuronSpikes(neuronSpikes>=start_times_ms(k) & neuronSpikes<=end_times_ms(k));
            neuronInd=max([round((neuronSpikesInWindow-start_times_ms(k))*ms2sample);ones(1,length(neuronSpikesInWindow))]);
            %make sure no bugs if spike is on the first Sample
            if ~isempty(neuronInd) && k==1 && neuronInd(1)==0 
                neuronInd(1)=1;
            end   
            binSpikes(i,(k-1)*window_samples+neuronInd)=true;
        end
    end
end

end

