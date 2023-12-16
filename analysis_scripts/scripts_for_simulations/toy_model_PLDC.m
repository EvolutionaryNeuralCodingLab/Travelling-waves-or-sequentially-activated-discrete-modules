%% Calculate dip p-value phase space (from simulation)
%(calculated in 210329)

%Notice: saveDirs not changed to relevant mats because different temporal
%settings have same name (but saved in different folder - temporalDir)


%clear all

%create spatial parameters
dxLim=[1 10];
sig_xLim=[1 10];
d_x=dxLim(1):0.1:dxLim(2);
sig_x=sig_xLim(1):0.1:sig_xLim(2);


plotLength=100;
nTimeSamples=25000; %this should be enough to avoid decretization errors in high spatial sigma,low dx situations


histogramEdges=0:(nTimeSamples/10):nTimeSamples;

% deltaT=1; %time different between times of peak height of the two guassians. Setting first peak to 0
% sigmaT=1; %gussians' std in time
temporalCombinations=[1 1 1 1 1; 0.5 0.75 1 2 3];
for temporalCombination=1:size(temporalCombinations,2)
    deltaT=temporalCombinations(1,temporalCombination);
    sigmaT=temporalCombinations(2,temporalCombination);
%     dip_p=zeros(length(sig_x),length(d_x));
    corrmat=zeros(length(sig_x),length(d_x));
    %         temporalDir=['/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/Phase space - DIP phase/deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) '/'];
    temporalDir=['/media/sil2/Data/Yuval O/Rotem Simulations/figures_for_paper/toy_model_PLDC/figures/Phase space - Correlation phase/deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) '/'];
    mkdir(temporalDir)

    t1=0;
    t2=t1+deltaT;
    t_avg=(t1+t2)/2;

    t0=t1;
    t_end=t2;
    for i=1:length(d_x)
        for j=1:length(sig_x)
            if any(i==1:20:length(d_x)) && any(j==1:20:length(sig_x))
                [deltaT,sigmaT,i,j]
                disp('printing')
                %                    saveFigs=1;
                saveFigs=0;
            else
                saveFigs=0;
            end
            deltaX=d_x(i);
            sigmaX=sig_x(j);
    %         deltaX=3;
    %         sigmaX=4;
            paramName=['deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) ' deltaX ' num2str(deltaX) ' sigmaX ' num2str(sigmaX)];
            x_avg=0;
            x1=-deltaX/2;
            x2=deltaX/2;

            %                 x=repmat(linspace(3*x1,3*x2,plotLength)',1,nTimeSamples);
            x=repmat(linspace(x1,x2,plotLength)',1,nTimeSamples);
            t=repmat(linspace(t0,t_end,nTimeSamples),plotLength,1); %total will be the time of the second peak plus two temporal stds

            v=exp(-(t-t1).^2/(2*sigmaT^2)-(x-x1).^2/(2*sigmaX^2))+exp(-(t-t2).^2/(2*sigmaT^2)-(x-x2).^2/(2*sigmaX^2));

            if saveFigs %plot signal at half time
                plot(x(:,round(nTimeSamples/2)),v(:,round(nTimeSamples/2))); %x is the same at all times, taking half time arbitrarily
                title([paramName ' Signal at half time'])
                xlabel('Spatial coordinate')
                ylabel('Signal')
                saveJpegAndFig(gcf,temporalDir,[paramName ' - Halftime Gaussians'],1);
                close gcf
            end

            [~,maximaTimes]=max(v,[],2);
            if saveFigs
                plot(maximaTimes,'.')
                title([paramName ' Maxima Latency'])
                xlabel('Horizontal Coordinate')
                ylabel('Time to Max')
                saveJpegAndFig(gcf,temporalDir,[paramName ' - Maxima Latency'],1);
                close gcf
            end
            R = corrcoef(x(:,1),maximaTimes);
            corrmat(j,i)=R(1,2);
            %                 [~, dip_p(i,j)] = hartigansdipsigniftest(sort(maximaTimes), 500);
            %in the phase diagram, dx (counted here as i) is columns
            %and sig_x (counted as j) is rows:
            [~, dip_p(j,i)] = hartigansdipsigniftest(sort(maximaTimes), 500);
            if saveFigs
                histogram(maximaTimes,histogramEdges)
                title(['Maxima Latencis Histogram - DIP pvalue ' num2str(dip_p(j,i))])
                saveJpegAndFig(gcf,temporalDir,[paramName ' - Histogram and p-value'],1);
                close gcf

                h=surf(v);
                set(h,'edgecolor','none');
                hold on
                scatter3(maximaTimes,1:plotLength,v(sub2ind(size(v),1:plotLength,maximaTimes')))
                saveJpegAndFig(gcf,temporalDir,[paramName ' - 3d wave'],1);
                close gcf
            end
        end
    end
    save([temporalDir 'dip values.mat'],'dip_p','d_x','sig_x','i','j','deltaT','sigmaT')
    save([temporalDir 'correlations.mat'],'corrmat','d_x','sig_x','i','j','deltaT','sigmaT')
%         load('/media/sil2/Literature/Projects/corplex/progress reports/meetings/next/Phase space - DIP phase/dip values.mat')

    %find numeric threshold
    [~,numThresh]=max(dip_p(:,1:end-1)>0.01 & dip_p(:,2:end)<=0.01,[],2);

    imagesc(d_x,sig_x,corrmat)
    hold on
    plot(d_x,d_x/2,'k','LineWidth',1)
    set(gca, 'YDir','normal')
    colorbar
    xlabel('\DeltaX [AU]')
    ylabel('\sigma_X [AU]')
    saveJpegAndFig(gcf,temporalDir,['deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) ' - correlation diagram - 0-1'],1);
%     saveJpegAndFig(gcf,temporalDir,['deltaT ' num2str(deltaT) ' sigmaT ' num2str(sigmaT) ' - correlation diagram'],1);
    close gcf

end
