function [comp_shared, comp_sharedIdx] = compareTopValues(medDist_top3idx, k_top3idx)
    % Compare top values from both datasets to find shared values
    comp_shared = zeros(size(medDist_top3idx));
    comp_sharedIdx = zeros(size(medDist_top3idx, 1), 1);
    
    for r = 1:size(medDist_top3idx, 1)
        tempMed = medDist_top3idx(r, :);
        tempK = k_top3idx(r, :);
        
        % Find shared values (zero where no match)
        comp_shared(r, :) = tempMed .* ismember(tempMed, tempK);
        
        % Find first match
        matchFound = false;
        for j = 1:length(tempMed)
            if ismember(tempMed(j), tempK) && tempMed(j) ~= 0
                comp_sharedIdx(r) = tempMed(j);
                matchFound = true;
                break;
            end
        end
        
        if ~matchFound
            comp_sharedIdx(r) = 0;
        end
    end
end