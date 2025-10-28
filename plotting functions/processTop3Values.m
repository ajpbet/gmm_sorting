function [k_top3idx, medDist_top3idx, comp_shared, comp_sharedIdx] = processTop3Values(kj_mat, polyID_M, medDist_cell, M_comp, exc_combPdf)
    % Process kj_mat to get top 3 values and their polyIDs
    [k_Mat, k_sortIdx, k_top3idx] = processKJMatrix(kj_mat, polyID_M);
    
    % Process medDist_cell with optional exclusion
    [medDist_sort, medDist_sortIdx, medDist_top3idx] = processMedDistMatrix(medDist_cell, polyID_M, M_comp, exc_combPdf);
    
    % Compare the top values from both sets
    [comp_shared, comp_sharedIdx] = compareTopValues(medDist_top3idx, k_top3idx);
end