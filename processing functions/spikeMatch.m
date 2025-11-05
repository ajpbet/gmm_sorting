function select_spike_match = spikeMatch(select_all, threshQ2, threshQ4, g, coeffs, folderSpikeMatch, basename)
    
    lSel = size(select_all,1);
    spike_idx = cell(lSel,1);
    
    for i = 1:lSel
        coeff_sel = select_all(i,3);
        gauss_select = select_all(i,4);
        mu = g{coeff_sel}.mu;
        sigma = sqrt(squeeze(g{coeff_sel}.Sigma(:)));
        mu_sel = mu(gauss_select);
        s_select = sigma(gauss_select);
        range_upp = mu_sel + 2*s_select;
        range_low = mu_sel - 2*s_select;
        spike_idx{i} = find(coeffs(:,coeff_sel) >= range_low & ...
                            coeffs(:,coeff_sel) <= range_upp);
    end
    
    % build pairwise overlap matrix
    n = numel(spike_idx);
    overlapMat = zeros(n);
    maxSizes = zeros(n);
    minSizes = zeros(n);
     for i = 1:n
        Ai = spike_idx{i};
        for j = i+1:n
            Aj = spike_idx{j};
            if isempty(Ai) || isempty(Aj)
                overlapMat(i,j) = 0;
            else
                inter = numel(intersect(Ai, Aj));
                szAi = numel(Ai);
                szAj = numel(Aj);
                
                minSz = min(szAj,szAi);
                minSizes(i,j) = minSz;
                maxSizes(i,j) = max(szAj,szAi);
                overlapMat(i,j) = inter / minSz;
            end
        end
    end
    
    % Calculate Normalized Minimum Distance to Boundary
    X = select_all(:,1);
    Y = select_all(:,2);
    
    x_min = min(X); x_max = max(X); x_range = x_max - x_min;
    y_min = min(Y); y_max = max(Y); y_range = y_max - y_min;
    
    if x_range == 0, x_range = 1; end
    if y_range == 0, y_range = 1; end
    
    X_norm = (X - x_min) / x_range;
    Y_norm = (Y - y_min) / y_range;
    
    m2_norm = threshQ2(1) * (x_range / y_range);
    b2_norm = (threshQ2(2) - y_min) / y_range;
    
    m4_norm = threshQ4(1) * (x_range / y_range);
    b4_norm = (threshQ4(2) - y_min) / y_range;
    
    dist_Q2 = abs(m2_norm * X_norm - Y_norm + b2_norm) ./ sqrt(m2_norm^2 + (-1)^2);
    dist_Q4 = abs(m4_norm * X_norm - Y_norm + b4_norm) ./ sqrt(m4_norm^2 + (-1)^2);
    
    x_inter = select_all(end,1);
    
    minDist_vec = zeros(n, 1);
    for k = 1:n
        if X(k) < x_inter
            minDist_vec(k) = dist_Q2(k);
        else
            minDist_vec(k) = dist_Q4(k);
        end
    end
    
    maxDistances = zeros(n);
    minDistances = zeros(n);
    
    for i = 1:n
        for j = i+1:n
            Di = minDist_vec(i);
            Dj = minDist_vec(j);
            
            minDistances(i,j) = min(Di, Dj);
            maxDistances(i,j) = max(Di, Dj);
        end
    end
    
    % Define the threshold for filtering
    threshold = 0.8;
    
    % PLOTTING DATA PREPARATION
    
    % Get the indices of the calculated upper triangle pairs
    upper_tri_indices = triu(ones(n), 1) == 1;
    
    % Extract the overlap values only for the calculated pairs
    overlap_vec = overlapMat(upper_tri_indices);
    
    % Identify the indices of pairs that meet the overlap threshold
    meets_threshold_indices = overlap_vec >= threshold;
    
    % Filter all data vectors based on the threshold
    X_minSizes = minSizes(upper_tri_indices);
    Y_maxSizes = maxSizes(upper_tri_indices);
    X_minDist = minDistances(upper_tri_indices);
    Y_maxDist = maxDistances(upper_tri_indices);
    
    X_minSizes_F = X_minSizes(meets_threshold_indices);
    Y_maxSizes_F = Y_maxSizes(meets_threshold_indices);
    X_minDist_F = X_minDist(meets_threshold_indices);
    Y_maxDist_F = Y_maxDist(meets_threshold_indices);

    % Filter the pair indices for labeling
    all_pair_indices = zeros(sum(upper_tri_indices(:)), 2); 
    k = 1;
    for i = 1:n
        for j = i+1:n
            all_pair_indices(k,:) = [i, j];
            k = k + 1;
        end
    end
    pair_indices_F = all_pair_indices(meets_threshold_indices, :);

    % Get the new total number of pairs to plot
    n_pairs_F = size(pair_indices_F, 1);
    
    % Generate colors based on the new, smaller number of pairs
    CoeffIDs = select_all(:, 3);
    GaussIDs = select_all(:, 4);
    
    all_hues = hsv(n_pairs_F);
    shuffled_order = randperm(n_pairs_F);
    pair_colors_F = all_hues(shuffled_order, :);
    
    % visualize the overlap
    hmap = figure;
    heatmap(overlapMat, 'ColorLimits', [0, 1], 'Title', sprintf('Pairwise Spike Overlap: %s', basename), ...
             'XLabel', 'Gauss Index', 'YLabel', 'Gauss Index');
    saveas(hmap, fullfile(folderSpikeMatch, ['overlap_' basename '.png']), 'png');
    close(hmap);
    % max and min spike count plot
    sizePlt = figure('Name', 'Spike Count Plot');
    ax_size = axes(sizePlt);
    hold(ax_size, 'on');
    
    for p = 1:n_pairs_F
        i_idx = pair_indices_F(p, 1);
        j_idx = pair_indices_F(p, 2);
        label = sprintf('C%dG%d, C%dG%d', CoeffIDs(i_idx), GaussIDs(i_idx), CoeffIDs(j_idx), GaussIDs(j_idx));

        % Scatter plot: Correct syntax uses 30 SizeData and color as positional args
        scatter(ax_size, X_minSizes_F(p), Y_maxSizes_F(p), 30, pair_colors_F(p,:), 'o', ...
                'filled', 'MarkerFaceAlpha', 0.5); 

        text(ax_size, X_minSizes_F(p), Y_maxSizes_F(p) + 0.03, label, ...
             'FontSize', 7, 'Color', pair_colors_F(p,:), 'HorizontalAlignment', 'center');
    end
    
    hold(ax_size, 'off');
    
    xlabel(ax_size, 'Minimum Spike Count minSz');
    ylabel(ax_size, 'Maximum Spike Count maxSz');
    title(ax_size, sprintf('Comparison of Spike Set Sizes for Overlapping Pairs: %s', basename));
    saveas(sizePlt, fullfile(folderSpikeMatch, ['size_' basename '.png']), 'png');
    close(sizePlt);

    % max and min boundary distance plot
    distPlt = figure('Name', 'Boundary Distance Plot');
    ax_dist = axes(distPlt);
    hold(ax_dist, 'on');
    
    for p = 1:n_pairs_F
        i_idx = pair_indices_F(p, 1);
        j_idx = pair_indices_F(p, 2);
        label = sprintf('C%dG%d, C%dG%d', CoeffIDs(i_idx), GaussIDs(i_idx), CoeffIDs(j_idx), GaussIDs(j_idx));
        
        % Scatter plot: Correct syntax uses 30 SizeData and color as positional args
        scatter(ax_dist, X_minDist_F(p), Y_maxDist_F(p), 30, pair_colors_F(p,:), 'o', ...
                'filled', 'MarkerFaceAlpha', 0.5); 
        
        text(ax_dist, X_minDist_F(p), Y_maxDist_F(p) + 0.03, label, ...
             'FontSize', 7, 'Color', pair_colors_F(p,:), 'HorizontalAlignment', 'center');
    end
    
    hold(ax_dist, 'off');
    
    xlabel(ax_dist, 'Minimum Normalized Distance to Boundary');
    ylabel(ax_dist, 'Maximum Normalized Distance to Boundary');
    title(ax_dist, sprintf('Pairwise Distance to Selection Boundary: %s', basename));
    saveas(distPlt, fullfile(folderSpikeMatch, ['dist_' basename '.png']), 'png'); 
    close(distPlt);

    % find groups of highly similar coefficients
    isSimilar = overlapMat >= threshold;
    G = digraph(isSimilar);
    groups = conncomp(G);
    
    % keep 1 representative per group
    numGroups = max(groups);
    reps = zeros(numGroups,1);
    for gnum = 1:numGroups
        idx = find(groups == gnum);
        lens = cellfun(@numel, spike_idx(idx));
        [~, bestIdx] = max(lens);
        reps(gnum) = idx(bestIdx);
    end
    select_spike_match = select_all(reps,:);
    
    numOriginal = size(select_all,1);
    numReduced  = size(select_spike_match,1);
    fprintf('Removed %d coefficients out of %d due to overlap\n', numOriginal - numReduced, numOriginal);
end