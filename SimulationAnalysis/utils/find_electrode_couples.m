function [two_scouts_results_sym, couples] = find_electrode_couples(EEGcap, two_scouts_results)
    n = size(EEGcap.Channel, 2);
    loc = zeros(n, 3);
    for i=1:n
        loc(i, :) = EEGcap.Channel(i).Loc;
    end

    mid_line_y_delta = 0.0007;
    loc_y = loc(:, 2);
    left_hemi_elec_b = loc_y > mid_line_y_delta;
    right_hemi_elec_b = loc_y < -mid_line_y_delta;
    fissure_elec_b = (~(left_hemi_elec_b | right_hemi_elec_b));
    
    % check all electrode are assigned only once
    electrode_assign = single(left_hemi_elec_b) + single(right_hemi_elec_b) + single(fissure_elec_b);
    if sum(electrode_assign == ones(n, 1)) ~= n
        warning("electrode assigninment problem")
    end
    
%     figure()
%     hold on
%     scatter_electrodes(loc(left_hemi_elec_b, :), 'b')
%     scatter_electrodes(loc(right_hemi_elec_b, :), 'r')
%     scatter_electrodes(loc(fissure_elec_b, :), 'k')
%     hold off


%     find couples
    couples = zeros(n, 1);
    for i=1:n
        if left_hemi_elec_b(i)
            ir = find__couple(loc, i);
            couples(i) = ir;

        elseif right_hemi_elec_b(i)
            il = find__couple(loc, i);
            couples(i) = il;

        else
            couples(i, :) = i;

        end
    
    end

%     figure()
%     scatter_electrodes_couples(loc, couples);

    two_scouts_results_sym = two_scouts_results;
    for rec=1: size(two_scouts_results, 1)
        elec_rec = two_scouts_results(rec, :).EEG_recordings{1};
        elec_rec_sym = zeros(size(elec_rec));
        for i=1:n
            elec_rec_sym(i, :) = elec_rec(i, :) + elec_rec(couples(i), :);
            
        end
        two_scouts_results_sym(rec, :).EEG_recordings{1} = elec_rec_sym;

    end

    
    
end

function I = find__couple(loc, i)
    loc_i = loc(i, :);
    loc_il_flip = [loc_i(1), -loc_i(2), loc_i(3)];
    dists = sqrt(sum((loc - loc_il_flip).^2, 2));
    [~, I] = min(dists);

end


function scatter_electrodes(loc, c)
    n = size(loc, 1);
    scatter3(loc(:, 1), loc(:, 2), loc(:, 3), c)

end


function scatter_electrodes_couples(loc, couples)

    n = size(loc, 1);
    scatter3(loc(:, 1), loc(:, 2), loc(:, 3))
    for i=1:n
        text(loc(i, 1), loc(i, 2), loc(i, 3), string(i) + "/" + string(couples(i)))
    end

end






