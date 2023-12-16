function [nChInWave,channels,times] = countContinousCrossings(allCrossings,En,maxTempDist,channels,times)
%COUNTCONTINOUSCROSSINGS is a recursive function that counts how many 
%channels have crossings that are connected (in space and time), i.e. how 
%many channels participate in a particular wave. 
%   Input:
%       - allCrossings: all channel's crossings times (one of the 4 
%         matrices given by getHilbertCrossings)
%       - En: channel layout 
%       - maxTempDist: the maximal time window allowed for two neighbors to
%       have crossings
%       - channels,times are arrays to keep track of the wave crossings 
%       found. AS INPUT, THE STARTING SEED SHOULD BE GIVEN (i.e. channels
%       should contain the seed channel and times the seed's time).
%   Output:
%       - nChInWave: Number of channels that had a countinous crossings 
%       from one another
%       - channels, times are arrays containing the channels that had
%       crossings and the time they had it

%TODO: - Maybe set  maxTempDist as default and allow change via varargin
%      - Currently looking only at horizontal and vertical neighbours.
%      Maybe diagonal should be added?

[currentPosI,currentPosJ]=find(En==channels(end));
currentCrossingTime=times(end);


neighborsSum=1; %There is always at least one - the seed

% check all neighbors. If they have crossing in temporal window, add 1 and 
% run function on them

for nextPosI=(currentPosI-1):(currentPosI+1)
    for nextPosJ=(currentPosJ-1):(currentPosJ+1)
        if (nextPosI==currentPosI && nextPosJ~=currentPosJ) || (nextPosI~=currentPosI && nextPosJ==currentPosJ) || (nextPosI==currentPosI && nextPosJ==currentPosJ) %Go only up,down,left,right, or stay in same channel
            if nextPosI>=1 && nextPosI<=size(En,1) && nextPosJ>=1 && nextPosJ<=size(En,2) && ~isnan(En(nextPosI,nextPosJ)) %make sure next channel is in layout and isn't a NaN
                nextChannel=En(nextPosI,nextPosJ);
                nextChannelCrossingsNoZeros=allCrossings(nextChannel,allCrossings(nextChannel,:)>0);
                closestCrossings=nextChannelCrossingsNoZeros(abs(nextChannelCrossingsNoZeros-currentCrossingTime)<maxTempDist);
                for i=1:numel(closestCrossings) %Find all closest crossings, which are closer than maxTempDist.
                    if all(~(channels==nextChannel & times==closestCrossings(i)))%next chanel hasn't apeared yet in the cluster
                        channels=[channels nextChannel];
                        times=[times closestCrossings(i)];
                        [nChInWaveNext,channels,times]=countContinousCrossings(allCrossings,En,maxTempDist,channels,times);
                        neighborsSum=neighborsSum+nChInWaveNext;
                    end
                end
            end
%         elseif  %address possibility that next crossing is in the same channel
%             nextChannel=channels(end); %same channel
%             nextChannelCrossingsNoZerosNoCurrent=allCrossings(nextChannel,allCrossings(nextChannel,:)>0 & allCrossings(nextChannel,:)~=currentCrossingTime);
%             %find the closest crossing, but not currentCrossingTime
%             closestCrossings=nextChannelCrossingsNoZerosNoCurrent(abs(nextChannelCrossingsNoZerosNoCurrent-currentCrossingTime)==min(abs(nextChannelCrossingsNoZerosNoCurrent-currentCrossingTime)));
%             for i=1:numel(closestCrossings) %There could be up to 2 closest crossing (one from each side). Need to address both becuase one might have appeared already in the cluster and the other not.
%                 if abs(closestCrossings(i)-currentCrossingTime)<maxTempDist && all(~(channels==nextChannel & times==closestCrossings(i)))%next chanel has a crossing within time window, and it hasn't apeared yet in the cluster
%                     channels=[channels nextChannel];
%                     times=[times closestCrossings(i)];
%                     [nChInWaveNext,channels,times]=countContinousCrossings(allCrossings,En,maxTempDist,channels,times);
%                     neighborsSum=neighborsSum+nChInWaveNext;
%                 end
%             end
        end
    end
end

nChInWave=neighborsSum;