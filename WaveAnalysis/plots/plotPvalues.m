function [f,PoD] = plotPvalues(p1,p2,varargin)
%PLOTPVALUES Summary of this function goes here
%   p1 is the left one, p2 is the right one
%   PoD is the probability of decrease (ratio of trials in which
%   p1>p2 out of all trials)

nSamples=length(p1);
if length(p2)~=nSamples
    error('p1 and p2 must be the same length')
end

plotLines=0;
plotBlueScatter=0;
KernelWidth=0.10;
leftColor='r';
rightColor='r';
diffColor ='k';
leftFaceAlpha=0.5;
rightFaceAlpha=0.5;
diffFaceAlpha=0.3;
scaling=0.3;
withmdn=1;
style=1;
p1Label='LFP';
p2Label='ALSA';
yLabelLeft='DIP Test P-Values';
yLabelRight='LFP-ALSA Pairwise Difference';
probOfDecrease=sum((p1-p2)>0)/nSamples;
Title='LFP and ALSA Comparison';
% xLabelText=['Prob of Decrease: ' num2str(probOfDecrease) ' (Out of ' num2str(nSamples) ' Trials)'];
plotMeanLine=1;
plotHalfOfLines=1;
valueLimits=1;
left_axis_color=[0 0 0];
right_axis_color=[0 0 0];
xlims=[0.5 2.5];

for i=1:2:length(varargin)
   eval([varargin{i} '=varargin{' num2str(i+1) '};']);
end


if plotHalfOfLines
    effectiveNSamples=max(1,round(nSamples/2));
    plotLineNums=randperm(nSamples,effectiveNSamples);
else
    effectiveNSamples=nSamples;
    plotLineNums=1:nSamples;
end


fig=figure;
set(fig,'defaultAxesColorOrder',[left_axis_color; right_axis_color]);


xlim(xlims)

if plotLines
    plot([ones(1,effectiveNSamples); 2*ones(1,effectiveNSamples)],[reshape(p1(plotLineNums),1,effectiveNSamples); reshape(p2(plotLineNums),1,effectiveNSamples)],'k');

    if plotBlueScatter
        hold on
        plot([ones(1,effectiveNSamples); 2*ones(1,effectiveNSamples)],[reshape(p1(plotLineNums),1,effectiveNSamples); reshape(p2(plotLineNums),1,effectiveNSamples)],'ob')
    end
    if plotMeanLine
        hold on
        plot([1 2],[mean(p1) mean(p2)],'b','LineWidth',5,'MarkerEdgeColor','r')
    end
    violinPositions=[1 2];
    xtickPositions=[1 2];
else
    violinPositions=[1.4 1.6];
    xtickPositions=[1.35 1.75];
%     plot(violinPositions,[mean(p1) mean(p2)],'b','LineWidth',3)
    p(1)=plot(violinPositions,[mean(p1) mean(p2)],'b','LineWidth',2);
    hold on
    p(2)=plot(violinPositions,[median(p1) median(p2)],'--k','LineWidth',2);
end

%plot violins
violin(violinPositions(1),p1, ...
       'KernelWidth', KernelWidth, ...
       'Rotation',    'vertical', ... 
       'side',        'left', ...
       'facecolor',   leftColor, ...
       'facealpha',   leftFaceAlpha, ...  
       'cutoff',      1e-3, ...
       'scaling',     scaling, ...
       'withmdn',    withmdn,...
       'style',         style,...
       'valueLimits',    valueLimits);

violin(violinPositions(2),p2, ...
       'KernelWidth', KernelWidth, ...
       'Rotation',    'vertical', ... 
       'side',        'right', ...
       'facecolor',   rightColor, ...
       'facealpha',   rightFaceAlpha, ...  
       'cutoff',      1e-3, ...
       'scaling',     scaling, ...
       'withmdn',    withmdn,...
       'style',         style,...
       'valueLimits',    valueLimits);

yyaxis right
ylim([-1 1])
yticks([-1 0 1])
ylabel(yLabelRight)

violin(xlims(2),(p1-p2), ...
       'KernelWidth', KernelWidth, ...
       'Rotation',    'vertical', ... 
       'side',        'left', ...
       'facecolor',   diffColor, ...
       'facealpha',   diffFaceAlpha, ...  
       'cutoff',      1e-3, ...
       'scaling',     scaling, ...
       'withmdn',    0,... %without median
       'style',         style,...
       'valueLimits',    valueLimits);
   
   
xticks(xtickPositions)
xticklabels({p1Label,p2Label})
set(gca,'TickLength',[0 0])

yyaxis left
ylim([0 1])
yticks([0 0.5 1])
ylabel(yLabelLeft)
title(Title)
% xlabel(xLabelText)
% annotation('textbox',[0.1 0 .1 .2],'String',boxText,'FitBoxToText','on','FontSize',8)

if ~plotLines 
    legend(p,'Mean','Median','Location','southwest')
end

if nargout
   f=fig;
end

if nargout>1
   PoD=probOfDecrease;
end

end

