function [] = saveJpegAndFig(f,savePath,imgName,seperateFigDir,varargin)
%SAVEJPEGANDFIG saves the figure f both as a matlab figure and a jpeg.
%   INPUT:
%       - f - handle of figure to save
%       - savePath - folder in which to save figure (must include a
%       filesep at the end
%       - imgName - name of the files (no extensions at the end)
%       - seperateFigDir - if true, saves the matlab fig file in a folder
%       withing savePath called 'figs' (also, creates it if it doesn't
%       exist). If false, saves both .fig and jpeg in savePath
%       - possible varargins (given as 'Key',value pairs)
%           - saveEPS - Also plots to EPS. Default is no. 

saveEPS=0;

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end


if seperateFigDir
%     figSaveDir=[savePath "figs" filesep];
    figSaveDir=savePath + "figs" + filesep;
    if ~exist(figSaveDir,'dir')
        mkdir(figSaveDir)
    end
else
    figSaveDir=savePath;
end

saveas(f,[savePath imgName '.jpg'])
if saveEPS
%     saveas(f,[savePath imgName],'epsc')
    print(f,'-r300',[savePath imgName '.eps'],'-depsc','-tiff')
end
% savefig(f,[figSaveDir imgName '.fig'])
savefig(f,figSaveDir + imgName + ".fig")
end

