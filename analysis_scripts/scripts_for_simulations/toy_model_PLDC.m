%% Calculate PLDC phase space 
%clear all

%create spatial parameters
dxLim=[1 10];
sig_xLim=[1 10];
d_x=dxLim(1):0.1:dxLim(2);
sig_x=sig_xLim(1):0.1:sig_xLim(2);

plotLength=100;
nTimeSamples=25000; %this should be enough to avoid decretization errors in high spatial sigma,low dx situations

deltaT=1; %time different between times of peak height of the two guassians. Setting first peak to 0
sigmaT=1; %gussians' std in time
corrmat=zeros(length(sig_x),length(d_x));

t1=0; % Time of first temporal peak
t2=t1+deltaT; % Time of second temporal peak
t_avg=(t1+t2)/2;

corrmat=zeros(length(sig_x),length(d_x));

t0=t1; % Simulation start time
t_end=t2; % Simulation end time

for i=1:length(d_x)
    for j=1:length(sig_x)
        deltaX=d_x(i);
        sigmaX=sig_x(j);
        x_avg=0;
        x1=-deltaX/2;
        x2=deltaX/2;

        x=repmat(linspace(x1,x2,plotLength)',1,nTimeSamples);
        t=repmat(linspace(t0,t_end,nTimeSamples),plotLength,1); %total will be the time of the second peak plus two temporal stds

        v=exp(-(t-t1).^2/(2*sigmaT^2)-(x-x1).^2/(2*sigmaX^2))+exp(-(t-t2).^2/(2*sigmaT^2)-(x-x2).^2/(2*sigmaX^2));

        [~,maximaTimes]=max(v,[],2);
        R = corrcoef(x(:,1),maximaTimes);
        corrmat(j,i)=R(1,2);
    end
end
figure
imagesc(d_x,sig_x,corrmat)
hold on
plot(d_x,d_x/2,'k','LineWidth',1)
set(gca, 'YDir','normal')
colorbar
xlabel('\DeltaX [AU]')
ylabel('\sigma_X [AU]')
title(['deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT)])



% set(gcf,'PaperPositionMode','auto');
% print 'PLDC_phase_space_dt2_s1' -dpdf -painters;