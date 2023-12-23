%% DIP-test p-values comparisons - Additional experiments summary

clear all
close all

%% Load final results for figure

% Load p-values comparison for ALSA
load([get_wave_analysis_code_base_path() 'precalculated_mats/DIP_values_additional_ALSA.mat'],...
    'pvalues','medAllDips')
pvalue_of = 'ALSA';

% % Load p-values comparison for first spikes
% load([get_wave_analysis_code_base_path() 'precalculated_mats/DIP_values_additional_ALSA.mat'],...
%     'pvalues','medAllDips')
% pvalue_of = 'spk';

%% Plot

figure;
hold on;
cellfun(@(x) plot([1 2],median(x,2),'LineWidth',2),pvalues,'UniformOutput',0);
ylim([0 1])
xlim([0.9 2.1])
xticks([1 2])
xticklabels({'LFP',pvalue_of})
plot([1 2],medAllDips,'--k','LineWidth',4);


set(gcf,'Units','centimeters','Position',[14.0000   12.9381    3.3    5.45]);
ylabel('DIP Test P-Values')
