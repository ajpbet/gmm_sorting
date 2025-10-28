function [medDist_cell_filtered, polyID_M_filtered] = excludeMcompValues(medDist_cell, polyID_M, M_comp)
    % Exclude values that are part of the combined PDF (M_comp)
    numRows = numel(M_comp);
    medDist_cell_filtered = cell(size(medDist_cell));
    polyID_M_filtered = cell(size(polyID_M));
    
    for r = 1:numRows
        vals = medDist_cell{r};
        ids = polyID_M{r};
        idx_exclude = M_comp{r};  % indices to exclude
        
        % Create mask to exclude specified indices and NaNs
        mask = ~ismember(ids, idx_exclude) & ~isnan(vals);
        
        medDist_cell_filtered{r} = vals(mask);
        polyID_M_filtered{r} = ids(mask);
    end
end