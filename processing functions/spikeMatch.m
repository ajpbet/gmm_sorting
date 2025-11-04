function select_spike_match = spikeMatch(select_all, g, coeffs)
    % This function removes redundant coefficients that select nearly
    % identical spikes (high overlap). It keeps one representative per group.
    %
    % select_all format: [medDist, kv, coeff, gauss]
    % g: Gaussian fit structs
    % coeffs: matrix of coefficient values

    lSel = size(select_all,1);
    spike_idx = cell(lSel,1);

    % --- Step 1: determine which spikes each Gaussian selects ---
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

    % --- Step 2: build pairwise overlap matrix (Jaccard similarity) ---
    n = numel(spike_idx);
    overlapMat = zeros(n);

    for i = 1:n
        Ai = spike_idx{i};
        for j = i+1:n
            Aj = spike_idx{j};

            if isempty(Ai) || isempty(Aj)
                overlapMat(i,j) = 0;
            else
                inter = numel(intersect(Ai, Aj));
                union_sz = numel(union(Ai, Aj));
                szAi = numel(Ai);
                szAj = numel(Aj);

                minSz = min(szAj,szAi);
                overlapMat(i,j) = inter / minSz;  % Jaccard similarity
            end
            overlapMat(j,i) = overlapMat(i,j);
        end
    end
    
    % visualize the overlap
    % colorbap plot of pariwise matrix of overlapMat
    % Visualize the overlap matrix using a heatmap
    figure;
    heatmap(overlapMat, 'ColorLimits', [0, 1], 'Title', 'Pairwise Spike Overlap', ...
             'XLabel', 'Spike Index', 'YLabel', 'Spike Index');
    % max and min distances to each boundary for each Ai v Aj for all pairs

    % compare max and min num of spikes for each pair

    

    % --- Step 3: find groups of highly similar coefficients ---
    threshold = 0.8;  % adjust as needed
    isSimilar = overlapMat >= threshold;
    G = graph(isSimilar);
    groups = conncomp(G);

    % --- Step 4: keep 1 representative per group (largest spike set) ---
    numGroups = max(groups);
    reps = zeros(numGroups,1);

    for gnum = 1:numGroups
        idx = find(groups == gnum);

        % choose the one that selects the *most* spikes
        lens = cellfun(@numel, spike_idx(idx));
        [~, bestIdx] = max(lens);
        reps(gnum) = idx(bestIdx);
    end

    select_spike_match = select_all(reps,:);
    
    % Print how many were removed
    numOriginal = size(select_all,1);
    numReduced  = size(select_spike_match,1);
    fprintf('Removed %d coefficients out of %d due to overlap\n', numOriginal - numReduced, numOriginal);

end
