function [k_select,k_lSel,kDist_vec] = KvKneeNoPkExc(k_sort, k_sortIdx,pg_init,g_init,xg_init)
% processMedDistKnee:
% Applies knee detection to each cell of medDist_sort using the same algorithm
% as processIdistKnee. Assumes each medDist_sort{k} is sorted ascending.

    numCells = numel(k_sort);
    k_select = cell(numCells, 1);
    k_lSel   = zeros(numCells, 1);

    % expands data into 2D vector, col 1 is medDist
    % and col 2 is coeff num (only used for plotting)
  % Calculate total number of elements for proper pre-allocation
    totalElements = sum(cellfun(@length, k_sort));
    k_vector = zeros(totalElements, 3);
    count = 0;
    
    % Expand all medDist to one vector
    for k = 1:numCells
        kV_k = k_sort{k};
        kv_Gidx = k_sortIdx{k};
        numElements = length(kV_k);
        
        [pk_vals, pk_locs, pk_widths] = findpeaks(pg_init{k}, xg_init{k}, 'WidthReference', 'halfheight');
        g = g_init{k};
        for i = 1:numElements
            mu = g.mu(kv_Gidx(i));
            std = g.Sigma(kv_Gidx(i));
            gauss_upp = mu+std;
            gauss_low = mu-std;
            inRangeIdx = pk_locs >= gauss_low & pk_locs <= gauss_upp;
            hasPeakInRange = any(inRangeIdx);
            
            if hasPeakInRange
                k_vector(count+1, 1) = kV_k(i);          % medDist values
                k_vector(count+1, 2) = k;                  % cell indices
                k_vector(count+1,3) = kv_Gidx(i);
                count = count + 1;
            end
        end        
    end
    
    % Sort by the first column (medDist values)
    [kVector_sort, kVectorIdx] = sort(k_vector(:, 1));
    kVector_sortCoeff = k_vector(kVectorIdx, 2);
    kVector_sortGauss = k_vector(kVectorIdx, 3);    
    
    kDist_vec = [kVector_sort,kVector_sortCoeff, ...
        kVector_sortGauss];
    
    sorted_vals = kVector_sort;
    sorted_idx  = kVector_sortCoeff;
    sorted_gauss = kVector_sortGauss;

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
    k_select{k}(:,1) = coeff_idx;           % coefficient indices
    k_select{k}(:,2) = kV_k(coeff_idx);% corresponding values
    k_select{k}(:,3) = gauss_idx;
    k_lSel(k) = length(coeff_idx);

end
