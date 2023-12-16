function channelData = convertMovieToChannels(movie,En)
%CONVERTMOVIETOCHANNELS converts a HeightXWidthXnFrames movie to ChxSamples
%matrix according to En channel map
%Usage: channelData = convertMovieToChannels(movie,En);
% FramesTot=size(movie,3);
nCh=max(En(:));
channelData=zeros(nCh,size(movie,3));
for i=1:nCh
    [row,col]=find(En==i);
    channelData(i,:)=movie(row,col,:);
end

end

