function [idist_kmatch,color_idistk,poly_match_idistK,ranking_Idist_match] = early_color_plot_proc(medDist_top3idx,k_top3idx)
 % first plot variables
    colors = cell(nRows,1);  % initialize 64x1, for example
    
    for i = 1:nRows
        matches_inTop3 = sum(ismember(medDist_top3idx(i,:), k_top3idx(i,:)));
        perfect_match = isequal(medDist_top3idx(i,:), k_top3idx(i,:));
    
        color = 'k';
        switch matches_inTop3
            case 0
                color = 'r';
            case 1
                color = 'y';
            case 2
                color = 'b';
            case 3
                if perfect_match
                    color = 'm';
                else
                    color = 'g';
                end
        end
        colors{i} = color;  % store it
    end
    
    % find which idist vals match a top 3 k val
     color_idistk = cell(nRows,1);
     idist_kmatch = zeros(nRows,1);
     ranking_Idist_match = zeros(nRows,1);
     poly_match_idistK = zeros(nRows,1);
     for i = 1:nRows
        vals = medDist_sort{i};       % numeric array
        ids  = medDist_sortIdx{i};    % numeric array
        
        nVals = numel(vals);
        for j = 1:nVals
            if ismember(ids(j), k_top3idx(i,:)) && ~isnan(vals(j))
                idist_kmatch(i) = vals(j);
                ranking_Idist_match(i) = j;
                poly_match_idistK(i) = ids(j);
                break;
            end
        end
     end
     for i = 1:length(ranking_Idist_match)
         color_idistk{i} = 'r'; % default for idist loc >3
         switch ranking_Idist_match(i)
             case 1
                color_idistk{i} = 'm';
             case 2
                color_idistk{i} = 'g';
             case 3
                color_idistk{i} = 'b';
         end
     end

end