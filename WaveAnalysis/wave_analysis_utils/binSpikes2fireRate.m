function spikingRate = binSpikes2fireRate(binSpikes,samplingFrequency,varargin)
%BINSPIKES2FIRERATE calculates the spiking rate (spike/s) according to
%binSpikes
%   Input:
%       -   binSpikes: logical nChannelsXnSamples with ones indicating
%       spiking events.
%       -   samplingFrequency: samples/sec of recording
%       -   Varargin (given as 'Name','Value' pairs):
%       	-   slidingWindowSize (1X1): size of the moving window in units
%       	of samples. Default is 10000 (0.5s for 20kHz)
%   Output:
%       - spikingRate: (nChXnSamples) the firing rate (spikes/sec)

slidingWindowSize=10000; %samples

for i=1:2:numel(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

normalization=slidingWindowSize/samplingFrequency;
avgFilt=ones(1,slidingWindowSize)/normalization;

spikingRate=conv2(binSpikes,avgFilt,'same');

end

