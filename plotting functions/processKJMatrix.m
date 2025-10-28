function [k_Mat, k_sortIdx, k_top3idx] = processKJMatrix(kj_mat, polyID_M)
    % Find maximum length for padding
    maxLen = max(cellfun(@numel, kj_mat));
    
    % Sort each numeric array in descending order
    k_Mat = cell(size(kj_mat));
    k_sortIdx = cell(size(kj_mat));
    
    for r = 1:numel(k_Mat)
        [vals_sorted, idx] = sort(kj_mat{r}, 'descend');  % sort original numeric values
        k_Mat{r} = vals_sorted;
        k_sortIdx{r} = polyID_M{r}(idx);        % map sorted indices to polyIDs
    end
    
    % Extract top 3 polyIDs for plotting/comparison
    k_top3idx = zeros(length(k_sortIdx), 3);
    for r = 1:length(k_sortIdx)
        nTop = min(3, length(k_sortIdx{r}));
        k_top3idx(r, 1:nTop) = k_sortIdx{r}(1:nTop);
    end
end