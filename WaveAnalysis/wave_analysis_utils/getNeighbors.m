function [nn,nnn] = getNeighbors(channel,En)
%GETNEIGHBORS returns arrays of 4 nearest and 4 next nearest neighboring
%channels according to En
%   Input:
%       - En: electrode layout
%   Output:
%       - nn: 4 (or less, if at edges) nearest neighbors
%       - nnn: 4 (or less, if at edges) of next nearest neighbors


[currentPosI,currentPosJ]=find(En==channel);

nn=[];
nnn=[];

for nextPosI=(currentPosI-1):(currentPosI+1)
    for nextPosJ=(currentPosJ-1):(currentPosJ+1)
        if nextPosI>=1 && nextPosI<=size(En,1) && nextPosJ>=1 && nextPosJ<=size(En,2) && ~isnan(En(nextPosI,nextPosJ)) %make sure next channel is in layout and isn't a NaN
%             if (nextPosI==currentPosI && nextPosJ~=currentPosJ) || (nextPosI~=currentPosI && nextPosJ==currentPosJ) || (nextPosI==currentPosI && nextPosJ==currentPosJ) %it is up, down, left or right
            if nextPosI~=currentPosI && nextPosJ~=currentPosJ %diagonal
                nnn=[nnn En(nextPosI,nextPosJ)];
            elseif ~(nextPosI==currentPosI && nextPosJ==currentPosJ) %up,down,left,right
                nn=[nn En(nextPosI,nextPosJ)];
            end
        end
    end
end

end

