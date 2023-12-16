function sampleCrossings=getCrossingsBySamples(crossings,hilbertAmps,varargin)
%getCrossingsBySamples returns the crossings of a trial in the form of
%nChannelsXnSamples matrix. The numeric value of each non-zero entry
%indicates the hilbert amplitude in the relevant channel at the time of 
%of the crossing. If hilbertAmp is an empty array, the value of each
%sample is 1. nSamples equals the sample of the last crossing in
%crossings.
%   Varargin:
%       - nSamples - if given, sampleCrossings will have nSamples columns,
%       filled with zeros after the last crossings (use this option to get
%       sampleCrossings to be the same size as other full trial variables
%       such as binSpikes given by getSpikeBinMatByChannel

nSamples=max(crossings(:));

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if isempty(hilbertAmps)
   hilbertAmps=crossings./crossings;
end

nCh=size(crossings,1);
sampleCrossings=zeros(nCh,nSamples);
for i=1:nCh
    nonZeroCrossings=crossings(i,:)>0;
    sampleCrossings(i,crossings(i,nonZeroCrossings))=hilbertAmps(i,nonZeroCrossings);
end

