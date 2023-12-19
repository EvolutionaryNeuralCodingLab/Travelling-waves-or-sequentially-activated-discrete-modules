function [f] = plot_phase_space(values_matrix, delta_Ts, delta_Xs, T_start_end,X_start_end,plot_significance,range,ax)
%PLOT_PHASE_SPACE Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('range','var') || isempty('range')
        range = [0 1];
    end
    if ~exist('ax','var')
        f = figure;
        ax = gca;
    else
        f = ancestor(ax,'figure');
    end
    imagesc(ax, values_matrix(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2)),range)
    colorbar()
    yticks(1:(T_start_end(2)-T_start_end(1)+1))
    yticklabels(num2cell(delta_Ts(T_start_end(1):T_start_end(2))))
    ylabel('Delta T / std(T)')
    xticks(1:(X_start_end(2)-X_start_end(1)+1))
    xticklabels(num2cell(delta_Xs(X_start_end(1):X_start_end(2))))
    xlabel('dipole dist / electrode dist')
    set(gca,'YDir','normal')
    if plot_significance
        hold( "on")
        [row,col] = ind2sub(size(values_matrix(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))),find(values_matrix(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))<0.1));
        scatter( col,row,'k','filled')
        [row,col] = ind2sub(size(values_matrix(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))),find(values_matrix(T_start_end(1):T_start_end(2),X_start_end(1):X_start_end(2))<0.05));
        scatter( col,row,'white')
    end
end

