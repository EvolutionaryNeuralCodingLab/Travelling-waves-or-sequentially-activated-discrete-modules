function simulatedPulses = simulateGaussians(layoutSize,gauss1Cov,gauss2Cov,pulseFrames,centerDist,tempOverlapInPulseFrames,varargin)
%SIMULATEGAUSSIANS creates a 3d representing an image series of two
%guassian pulses
%
%   Input:
%
%   Channel layout is layoutSizeXlayoutSize
%   gauss1Cov,gauss2Cov are the 1X1 or 2X2 are the covariance matrices of the first and second Guassians. 
%   If gaussCov is 1x1 then matrix is this value times the unit matrix
%   pulseFrames are the total frame number from the apearance of a gauss
%   to its disappearance.
%   centerDist is [distInX distInY], the distance in x,y from the center
%   of the first gaussian to the center of the next. Default units are
%   channels/pixels, but can be changed using 'distUnits' varargin.
%   tempOverlapInPulseFrames is a the percentage of temporal overlap between
%   the two pulses (i.e. takes a number between 0 and 1)
%   Possible varargins pairs ('Name',value):
%       - temporal and spatial positions of first pulse,x1,y1,t1
%       - temporalFunc is the time function multiplying the gaussian.
%         Should be 1xpulseFrames long. Default is 
%         sin(linspace(0,pi,pulseFrames)).^2
%       - distUnits: the units in which centerDist is given. Defult is
%       'Channels' (pixels) but can be changes to 'Sigma' if gauss1Cov is
%       a 1x1 scalar, in which case distance will be in units of
%       sqrt(gauss1Cov) (i.e. the std of the first gaussian)
%       - saveSpatialShapes: logical. save the simulated gaussian shapes
%       (1,2 and the sum of them) if true. defult is false. if true,
%       saveDir should be given
%       - saveDir - folder in which to save gaussian shapes. required for when 
%       saveSpatialShapes is true. Must have filesep at the end
%   Output:
%   
%   simulatedPulses is an layoutSizeXlayoutSizeXt_tot (default t_tot is 
%   the time the second pulse ends). 


%Define variables
x1=layoutSize/4;
y1=layoutSize/4;
t1=1;

distUnits='Channels';

temporalFunc=sin(linspace(0,pi,pulseFrames)).^2;

saveSpatialShapes=0; 

%overwrite any of these using varargin pairs

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end

if saveSpatialShapes && ~exist('saveDir','var')
    error('For saving spatial shapes (saveSpatialShapes=1) saveDir must be also given')
end

if strcmp(distUnits,'Sigma')
    if ~numel(gauss1Cov)==1
        error('For ditsUnits to be ''Sigma'', gauss1Cov must be 1x1')
    end
    centerDist=centerDist*sqrt(gauss1Cov);
elseif ~strcmp(distUnits,'Channels')
    error('Wrong input for ''distUnits''! Must be either ''Sigma'' or ''Channels''')
end

if det(gauss1Cov)==0 || det(gauss2Cov)==0
    error('A covariance matrix is singular! Enter non singular covariance matrices only')
end
if numel(gauss1Cov)==1
    gauss1Cov=eye(2)*gauss1Cov;
end
if numel(gauss2Cov)==1
    gauss2Cov=eye(2)*gauss2Cov;
end

x2=x1+centerDist(1);
y2=y1+centerDist(2);
t2=t1+pulseFrames-pulseFrames*tempOverlapInPulseFrames;
t_tot=t2+pulseFrames-1;

mu1=[x1,y1];
mu2=[x2,y2];

simulatedPulses=zeros(layoutSize,layoutSize,t_tot);
[X,Y]=meshgrid(1:layoutSize,1:layoutSize);
gauss1=zeros(layoutSize);
gauss2=zeros(layoutSize);

for i=1:layoutSize
   for j=1:layoutSize
        pos=[X(i,j),Y(i,j)];
        gauss1(i,j)=exp(-(pos-mu1)*(gauss1Cov^-1)*(pos-mu1)'/2);
        gauss2(i,j)=exp(-(pos-mu2)*(gauss2Cov^-1)*(pos-mu2)'/2);
   end
end

if saveSpatialShapes
    imagesc(gauss1)
    saveJpegAndFig(gcf,saveDir,'Gauss 1',1);
    close gcf
    imagesc(gauss2)
    saveJpegAndFig(gcf,saveDir,'Gauss 2',1);
    close gcf
    imagesc(gauss1+gauss2)
    saveJpegAndFig(gcf,saveDir,'Gauss 1+2',1);
    close gcf
% else
%     figure
%     imagesc(gauss1)
%     figure
%     imagesc(gauss2)
end

pulse1=repmat(gauss1,1,1,pulseFrames).*shiftdim(repmat(temporalFunc,layoutSize,1,layoutSize),2);
pulse2=repmat(gauss2,1,1,pulseFrames).*shiftdim(repmat(temporalFunc.^2,layoutSize,1,layoutSize),2);

simulatedPulses(:,:,t1:(t1+pulseFrames-1))=pulse1;
simulatedPulses(:,:,t2:(t2+pulseFrames-1))=simulatedPulses(:,:,t2:(t2+pulseFrames-1))+pulse2;


end

