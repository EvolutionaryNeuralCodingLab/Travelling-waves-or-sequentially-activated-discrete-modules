%clear all

layoutSize=12;
En=reshape(1:(layoutSize^2),layoutSize,layoutSize);

% radial waves

pulseFrames=100;
distInSigmas=[0 0];

tempOverlapInPulseFrames=0.7;
gaussSigma=[3 6];
radialWave=simulateGaussians(layoutSize,gaussSigma(1)^2,gaussSigma(2)^2,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'x1',layoutSize/2,'y1',layoutSize/2,'distUnits','Sigma');

HT=hilbert(squeeze(convertMovieToChannels(radialWave,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(radialWave,3)]; 
plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),[],'plotColorBar',0)
radial_PLDC = calcDistanceAndPhaseLatencyCorrelation(1:numel(En),crossings{1}(:,1),En);
set(gcf,'Position',[ 23.7331   13.4938    2.4342    3.4660]);

plotFilmstripSimulated(radialWave,[ceil(pulseFrames/4),pulseFrames],En,crossings{1},hilbertAmps{1});


% two gaussians slightly overlapping (time+space)
layoutSize=12;
gaussSigma=6;
pulseFrames=100;
distInSigmas=[1 1];
tempOverlapInPulseFrames=0.7;
waveData=simulateGaussians(layoutSize,gaussSigma^2,gaussSigma^2,pulseFrames,distInSigmas,tempOverlapInPulseFrames,'distUnits','Sigma');
HT=hilbert(squeeze(convertMovieToChannels(waveData,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(waveData,3)]; 

plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),[],'plotColorBar',0)

plotFilmstripSimulated(waveData,[ceil(pulseFrames/4),pulseFrames],En,crossings{1},hilbertAmps{1});


% spiral waves 
pulseFrames=100;
distBetweenCenters=[0 4];
guassianAngles=[-pi/4 +pi/4];
rotate1=[cos(guassianAngles(1)) -sin(guassianAngles(1));sin(guassianAngles(1)) cos(guassianAngles(1))];
rotate2=[cos(guassianAngles(2)) -sin(guassianAngles(2));sin(guassianAngles(2)) cos(guassianAngles(2))];

verticalEllips=[8 0; 0 15];
cov1=rotate1*verticalEllips*rotate1^-1;
cov2=rotate2*verticalEllips*rotate2^-1;
tempOverlapInPulseFrames=0.7;
gauss1centerPosition=[layoutSize/2+2,layoutSize/3];
spiralWave=simulateGaussians(layoutSize,cov1,cov2,pulseFrames,distBetweenCenters,tempOverlapInPulseFrames,'x1',gauss1centerPosition(1),'y1',gauss1centerPosition(2));

HT=hilbert(squeeze(convertMovieToChannels(spiralWave,En))').';
HTabs=abs(HT);
HTangle=angle(HT);

[crossings,hilbertAmps] = getHilbertCrossings(HTabs,HTangle);

startEndWave=[1 size(spiralWave,3)]; 

plotCrossingsPhysical(crossings{1},startEndWave,flipud(En),[],'plotColorBar',0)

plotFilmstripSimulated(spiralWave,[ceil(pulseFrames/4),pulseFrames],En,crossings{1},hilbertAmps{1});
