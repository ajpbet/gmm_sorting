function [coeff_idist,idist_select] = idist_knee_select(idist_kmatch)
    [sorted_idist, ind_idist] = sort(idist_kmatch);  

    % --- Keep only top maxIdist candidates for knee detection ---
    minIdist = 1;
    maxIdist = 64;
    top_candidates = sorted_idist(end-maxIdist+1:end);
    top_indices    = ind_idist(end-maxIdist+1:end);
    
    ncoeff = length(sorted_idist);
    maxA_idist = max(sorted_idist);
    
    % --- Sliding-window slope computation ---
    nd = 10;
    if ncoeff >= nd
        d_idist = (top_candidates(nd:end) - top_candidates(1:end-nd+1)) ...
                   / maxA_idist * ncoeff / nd;
        all_above1_idist = find(d_idist >= 1);
    else
        all_above1_idist = [];
    end
    
    % --- Knee detection ---
    if numel(all_above1_idist) >= 2
        aux2_idist = diff(all_above1_idist);
        temp_bla_idist = conv(aux2_idist(:), [1 1 1]/3);   % smooth differences
        temp_bla_idist = temp_bla_idist(2:end-1);
        temp_bla_idist(1) = aux2_idist(1);
        temp_bla_idist(end) = aux2_idist(end);
    
        % First position where 3 consecutive differences == 1
        thr_knee_diff_idist = all_above1_idist(find(temp_bla_idist(2:end) == 1, 1)) + nd/2;
    
        % Number of inputs to select
        inputs_idist = maxIdist - thr_knee_diff_idist + 1;
    else
        % Fallback if no steep slope detected
        inputs_idist = minIdist;
    end
    
    % --- Select top coefficients ---
    coeff_idist = top_indices(end:-1:end-inputs_idist+1);   % descending order
    idist_select(:,1) = coeff_idist;                           % coefficient indices
    idist_select(:,2) = idist_kmatch(coeff_idist);            % corresponding IDist values
    

end