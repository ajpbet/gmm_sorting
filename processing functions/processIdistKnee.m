function [idist_select, inputs_idist, idist_kmatch_lSel, sorted_idist, ind_idist] = ...
         processIdistKnee(idist_kmatch)
% PROCESSIDISTKNEE: Selects top IDist coefficients using a knee-detection method
%
% INPUTS:
%   idist_kmatch  - Vector of IDist values (e.g., per coefficient)
%   minIdist      - Minimum number of coefficients to select (default = 1)
%   maxIdist      - Maximum number of coefficients considered for knee detection (default = 30)
%   nd            - Window size for slope detection (default = 10)
%
% OUTPUTS:
%   idist_select        - Nx2 matrix: [coeff_index, IDist_value] of selected coefficients
%   inputs_idist        - Number of coefficients selected after knee detection
%   idist_kmatch_lSel   - Count of selected coefficients (redundant but useful)
%   sorted_idist        - Sorted IDist values (ascending)
%   ind_idist           - Original indices corresponding to sorted IDist
%
% EXAMPLE:
%   [idist_select, inputs_idist] = processIdistKnee(idist_kmatch);
%
%   This selects the top coefficients based on knee detection in the
%   sorted IDist curve and returns their indices and values.

    % --- Handle defaults ---
    nd = 10;
    maxIdist = 64;
    minIdist = 1;

    % --- Sort IDist values ---
    [sorted_idist, ind_idist] = sort(idist_kmatch);

    non_zero_mask = sorted_idist ~= 0;
    sorted_idist = sorted_idist(non_zero_mask);
    ind_idist = ind_idist(non_zero_mask);
    ncoeff = length(sorted_idist);
    maxA_idist = max(sorted_idist);

    % --- Extract top candidates ---
    if ncoeff < maxIdist
        maxIdist = ncoeff;
    end
    top_candidates = sorted_idist(end - maxIdist + 1:end);
    top_indices    = ind_idist(end - maxIdist + 1:end);

    % --- Sliding-window slope computation ---
    if ncoeff >= nd
        d_idist = (top_candidates(nd:end) - top_candidates(1:end-nd+1)) ...
                    / maxA_idist * ncoeff / nd;
        all_above1_idist = find(d_idist >= 1);
    else
        all_above1_idist = [];
    end

    % --- Knee detection logic ---
    if numel(all_above1_idist) >= 2
        aux2_idist = diff(all_above1_idist);
        temp_bla_idist = conv(aux2_idist(:), [1 1 1]/3);   % smooth differences
        temp_bla_idist = temp_bla_idist(2:end-1);
        temp_bla_idist(1) = aux2_idist(1);
        temp_bla_idist(end) = aux2_idist(end);

        % First position where 3 consecutive diffs == 1
        thr_knee_diff_idist = all_above1_idist(find(temp_bla_idist(2:end) == 1, 1)) + nd/2;

        % Number of coefficients to select
        if isempty(thr_knee_diff_idist)
            inputs_idist = minIdist;
        else
            inputs_idist = maxIdist - thr_knee_diff_idist + 1;
        end
    else
        % Fallback if no steep slope detected
        inputs_idist = minIdist;
    end

    % --- Clamp selection count ---
    inputs_idist = max(min(inputs_idist, maxIdist), minIdist);

    % --- Select top coefficients (descending order) ---
    coeff_idist = top_indices(end:-1:end-inputs_idist+1);
    idist_select = [coeff_idist(:), idist_kmatch(coeff_idist(:))];

    % --- Output number of selected coefficients ---
    idist_kmatch_lSel = length(idist_select(:,1));

end
