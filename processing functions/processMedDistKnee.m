function [medDist_select,medDist_lSel,medDist_vec] = processMedDistKnee(medDist_sort, medDist_sortIdx)
% processMedDistKnee:
% Applies knee detection to each cell of medDist_sort using the same algorithm
% as processIdistKnee. Assumes each medDist_sort{k} is sorted ascending.

    numCells = numel(medDist_sort);
    medDist_select = cell(numCells, 1);
    medDist_lSel   = zeros(numCells, 1);

    % expands data into 2D vector, col 1 is medDist
    % and col 2 is coeff num (only used for plotting)
  % Calculate total number of elements for proper pre-allocation
    totalElements = sum(cellfun(@length, medDist_sort));
    medDist_vector = zeros(totalElements, 3);
    count = 0;
    
    % Expand all medDist to one vector
    for k = 1:numCells
        medDist_k = medDist_sort{k};
        medDist_Gidx = medDist_sortIdx{k};
        numElements = length(medDist_k);
        
        % Fill the appropriate section
        indices = count + (1:numElements);
        medDist_vector(indices, 1) = medDist_k;          % medDist values
        medDist_vector(indices, 2) = k;                  % cell indices
        medDist_vector(indices,3) = medDist_Gidx;
        count = count + numElements;
    end
    
    % Sort by the first column (medDist values)
    [medVector_sort, medVectorIdx] = sort(medDist_vector(:, 1));
    medVector_sortCoeff = medDist_vector(medVectorIdx, 2);
    medVector_sortGauss = medDist_vector(medVectorIdx, 3);    
    
    medDist_vec = [medVector_sort,medVector_sortCoeff, ...
        medVector_sortGauss];
    
    sorted_vals = medVector_sort;
    sorted_idx  = medVector_sortCoeff;
    sorted_gauss = medVector_sortGauss;

    % --- Keep only top max candidates for knee detection ---
    minVal = 1;
    maxVal = 64;
    top_candidates = sorted_vals(end-maxVal+1:end);
    top_indices    = sorted_idx(end-maxVal+1:end);
    top_gauss      = sorted_gauss(end-maxVal+1:end);
    ncoeff = length(sorted_vals);
    maxA = max(sorted_vals);
    
    % --- Sliding-window slope computation ---
    nd = 10;
    if ncoeff >= nd
        d_vals = (top_candidates(nd:end) - top_candidates(1:end-nd+1)) ...
                 / maxA * ncoeff / nd;
        all_above1 = find(d_vals >= 1);
    else
        all_above1 = [];
    end
    
    % --- Knee detection ---
    if numel(all_above1) >= 2
        aux2 = diff(all_above1);
        temp_bla = conv(aux2(:), [1 1 1]/3); % smooth differences
        temp_bla = temp_bla(2:end-1);
        temp_bla(1) = aux2(1);
        temp_bla(end) = aux2(end);
    
        % First position where 3 consecutive differences == 1
        thr_knee_diff = all_above1(find(temp_bla(2:end) == 1, 1)) + nd/2;
    
        % Number of inputs to select
        inputs = maxVal - thr_knee_diff + 1;
    else
        % Fallback if no steep slope detected
        inputs = minVal;
    end
    
    % --- Select top coefficients ---
    coeff_idx = top_indices(end:-1:end-inputs+1);
    gauss_idx = top_gauss(end:-1:end-inputs+1);
    medDist_select{k}(:,1) = coeff_idx;           % coefficient indices
    medDist_select{k}(:,2) = medDist_k(coeff_idx);% corresponding values
    medDist_select{k}(:,3) = gauss_idx;
    medDist_lSel(k) = length(coeff_idx);

end
