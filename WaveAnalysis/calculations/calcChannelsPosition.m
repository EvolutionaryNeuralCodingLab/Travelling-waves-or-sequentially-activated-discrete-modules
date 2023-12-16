function chPos = calcChannelsPosition(En)
%CALCCHANNELSPOSITION calculates the position of channels according to the
%electrode layout En
%   INPUT:
%       En - Electrode Layout
%   OUTPUT:
%   chPos - nChx2 of all channels' positions. chPos(i,1) and chPos (i,2)
%   are the row and column of channel i respectively.
nCh=max(En(:));

chPos=zeros(nCh,2);
for i=1:nCh
    [chPos(i,1),chPos(i,2)]=find(En==i);
end
end

