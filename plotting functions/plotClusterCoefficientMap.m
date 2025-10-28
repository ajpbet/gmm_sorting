function plotClusterCoefficientMap(g, medDist_sort, medDist_sortIdx, wd_coeff, spikes, cluster_times, coeff_vals, channelNum, folderSpike)
    % PLOTCLUSTERCOEFFICIENTMAP:
    % Visualize per-cluster percentage of spikes captured by top Gaussians in each coefficient.
    % Uses medDist_sort and medDist_sortIdx to select top Gaussians.
    %
    % Inputs:
    %   g               - cell array of GMM structs (g{coeff_num}.mu, Sigma)
    %   medDist_sort    - cell array, sorted metric values for Gaussian quality (filtered already)
    %   medDist_sortIdx - cell array, corresponding Gaussian indices (may include zeros)
    %   wd_coeff        - [spikes x coeffs] waveform-derived coefficients
    %   spikes          - actual spike waveforms (spike count x features)
    %   cluster_times   - [spike x 1] cluster ID assignments
    %   coeff_vals      - array of coefficient numbers (e.g., 1:8)
    %   channelNum      - string/number for labeling
    %   folderSpike     - output folder for saving results

    fprintf('--- Running Cluster-Coefficient Map ---\n');

    all_clusters = unique(cluster_times(:,1));
    num_clusters = numel(all_clusters);
    num_coeffs   = numel(coeff_vals);

    % Total spikes per cluster for normalization
    total_spikes_per_cluster = zeros(num_clusters,1);
    for c = 1:num_clusters
        total_spikes_per_cluster(c) = sum(cluster_times(:,1) == all_clusters(c));
    end

    % Initialize matrices
    percentage_matrix = zeros(num_clusters, num_coeffs);

    % -------------------------
    % Local behavior flag:
    %   true  -> use only the highest (first) Gaussian in medDist_sortIdx
    %   false -> use top 3 (or fewer if not available)
    use_first_only = true;
    % -------------------------

    % === MAIN LOOP: each coefficient ===
    for coeff_idx = 1:num_coeffs
        coeff_num = coeff_vals(coeff_idx);

        % safety: ensure g has this coeff
        if coeff_num > numel(g) || isempty(g{coeff_num}) || isempty(g{coeff_num}.mu)
            fprintf('Coeff %d -> no GMM data, skipping.\n', coeff_num);
            continue;
        end

        gmm = g{coeff_num};

        % Retrieve medDist entries using coeff_num as index into cells
        if coeff_num <= numel(medDist_sort) && ~isempty(medDist_sort(coeff_num)) ...
                && coeff_num <= numel(medDist_sortIdx) && ~isempty(medDist_sortIdx(coeff_num))
            if use_first_only
                % take only the first (highest) entry if it exists
                candidate = medDist_sortIdx(coeff_num,1);
                top_idx = candidate;
            else
                % top-3 (or fewer if not enough)
                nWanted = 3;
                avail = numel(medDist_sortIdx(coeff_num));
                nUse = min(nWanted, avail);
                top_idx = medDist_sortIdx(coeff_num,1:nUse);
            end
            % Ensure top_idx is a row vector and numeric
            top_idx = top_idx(:).';
            % Remove zeros and OOB indices (must map to indices into gmm.mu)
            top_idx = top_idx(top_idx > 0 & top_idx <= numel(gmm.mu));
        else
            top_idx = [];
        end

        if isempty(top_idx)
            fprintf('Coeff %d → No valid Gaussians found, skipping.\n', coeff_num);
            continue;
        end

        fprintf('Coeff %d using Gaussians: %s\n', coeff_num, mat2str(top_idx));

        % Gather spikes that fall into any of the selected Gaussians' ranges
        spike_mask = false(size(wd_coeff,1),1);
        for g_idx = top_idx
            % defensive: ensure g_idx is within bounds
            if g_idx < 1 || g_idx > numel(gmm.mu)
                continue;
            end

            mu_g    = gmm.mu(g_idx);
            % Sigma may be stored as 1x1xK or DxDxK; handle both
            S = gmm.Sigma;
            if ndims(S) == 3
                sigma_g = sqrt(squeeze(S(:,:,g_idx)));
            else
                sigma_g = sqrt(squeeze(S(g_idx,:)));
            end

            % if sigma_g is vector, loop elements; otherwise treat scalar
            if numel(mu_g) == numel(sigma_g)
                % handle element-wise ranges
                in_range = false(size(wd_coeff,1),1);
                for elem = 1:numel(mu_g)
                    range_st  = mu_g(elem) - 2*sigma_g(elem);
                    range_end = mu_g(elem) + 2*sigma_g(elem);
                    spike_vals = wd_coeff(:, coeff_num);
                    spike_vals = spike_vals(:);
                    in_range = in_range | (spike_vals >= range_st & spike_vals <= range_end);
                end
            else
                % fallback (treat mu_g/sigma_g as scalars)
                range_st  = mu_g - 2*abs(sigma_g);
                range_end = mu_g + 2*abs(sigma_g);
                spike_vals = wd_coeff(:, coeff_num);
                spike_vals = spike_vals(:);
                in_range = (spike_vals >= range_st & spike_vals <= range_end);
            end

            spike_mask = spike_mask | in_range;  % union of all selected Gaussians
        end

        spikeIdx = find(spike_mask);
        if isempty(spikeIdx)
            fprintf('Coeff %d → no spikes found in selected Gaussian ranges.\n', coeff_num);
            continue;
        end

        spike_clusters = cluster_times(spikeIdx,1);

        % Compute per-cluster percentage for this coefficient
        for c = 1:num_clusters
            cluster_id = all_clusters(c);
            num_in_cluster = sum(spike_clusters == cluster_id);
            if total_spikes_per_cluster(c) > 0
                percentage_matrix(c, coeff_idx) = (num_in_cluster / total_spikes_per_cluster(c)) * 100;
            else
                percentage_matrix(c, coeff_idx) = 0;
            end
        end
    end

    % === Visualization ===
    fig = figure('Name','Cluster-Coefficient Percentage Heatmap', ...
        'NumberTitle','off','Position',[100 100 1400 900]);
    set(fig, 'WindowStyle', 'docked');

    ax = axes('Parent', fig);
    imagesc(percentage_matrix, [0 100]);
    axis tight;
    axis equal;

    % Red→Green colormap (low=red, high=green)
    cmap = flipud(turbo);
    colormap(ax, cmap);
    colorbar;
    caxis([0 100]);

    % Labels
    xlabel('Coefficient');
    ylabel('Cluster');
    xticks(1:num_coeffs);
    yticks(1:num_clusters);
    xticklabels(arrayfun(@(x) sprintf('C%d', x), coeff_vals, 'UniformOutput', false));
    yticklabels(arrayfun(@(x) sprintf('Clust %d', x), all_clusters, 'UniformOutput', false));
    title(sprintf('Cluster Spikes Captured per Coefficient (Channel %s)', num2str(channelNum)), 'FontWeight','bold');

    % Annotate percentage values
    for c = 1:num_clusters
        for coeff_idx = 1:num_coeffs
            val = percentage_matrix(c, coeff_idx);
            if val > 0
                if val > 50
                    text_color = 'white';
                else
                    text_color = 'black';
                end
                text(coeff_idx, c, sprintf('%.1f%%', val), ...
                     'HorizontalAlignment','center', ...
                     'VerticalAlignment','middle', ...
                     'FontSize',8,'FontWeight','bold','Color',text_color);
            end
        end
    end

    % Save output
    if ~exist(folderSpike, 'dir')
        mkdir(folderSpike);
    end
    filename = fullfile(folderSpike, sprintf('ch%s_cluster_coefficient_heatmap.png', num2str(channelNum)));
    exportgraphics(fig, filename, 'Resolution', 300);

    fprintf('Saved cluster-coefficient heatmap to %s\n', filename);
end
