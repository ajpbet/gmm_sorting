function [medDist_sort, medDist_sortIdx, medDist_top3idx] = processMedDistMatrix(medDist_cell, polyID_M, M_comp, exc_combPdf)
    % Apply exclusion first if requested
    if exc_combPdf
        [medDist_cell, polyID_M] = excludeMcompValues(medDist_cell, polyID_M, M_comp);
    end
    
    % Sort medDist values
    medDist_sort = cell(size(medDist_cell));
    medDist_sortIdx = cell(size(medDist_cell));
    
    for r = 1:numel(medDist_cell)
        vals = medDist_cell{r};
        ids = polyID_M{r};
        
        % Exclude NaNs
        mask = ~isnan(vals);
        vals_valid = vals(mask);
        ids_valid = ids(mask);
        
        % Sort the valid values
        [vals_sorted, idx] = sort(vals_valid, 'descend');
        
        medDist_sort{r} = vals_sorted;          % store sorted medDist values
        medDist_sortIdx{r} = ids_valid(idx);    % map sorted indices to polyIDs
    end
    
    % Extract top 3 polyIDs
    medDist_top3idx = zeros(numel(medDist_sort), 3);
    for r = 1:numel(medDist_sort)
        vals = medDist_sort{r};
        ids = medDist_sortIdx{r};
        nTop = min(3, numel(vals));
        
        if nTop > 0
            medDist_top3idx(r, 1:nTop) = ids(1:nTop);
        end
    end
end