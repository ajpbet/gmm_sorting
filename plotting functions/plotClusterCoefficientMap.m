function plotClusterCoefficientMap(g, polyID_M, wd_coeff, spikes, cluster_times, coeff_vals, channelNum, folderSpike)
    % PLOTCLUSTERCOEFFICIENTMAP:
    %   Visualize per-cluster Gaussian performance per coefficient.
    %   Figure 1: max percentage over all Gaussians (color = %)
    %   Figure 2: Gaussian index that achieved that max (color = G#)
    %
    % Inputs:
    %   g           - cell array of GMM structs
    %   polyID_M    - cell array, each cell contains the Gaussian indices present for that coeff
    %   wd_coeff    - [spikes x coeffs]
    %   spikes      - [spikes x features]
    %   cluster_times - [spike x 1] cluster IDs
    %   coeff_vals  - array of coefficients (1:8)
    %   channelNum  - label or number
    %   folderSpike - output folder path

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
    winner_matrix     = nan(num_clusters, num_coeffs);  % NaN for no spikes / no relevant Gaussians

    % --- MAIN LOOP over coefficients ---
    for coeff_idx = 1:num_coeffs
        coeff_num = coeff_vals(coeff_idx);

        % Skip empty GMMs
        if coeff_num > numel(g) || isempty(g{coeff_num}) || isempty(g{coeff_num}.mu)
            fprintf('Coeff %d → no GMM data.\n', coeff_num);
            continue;
        end
        gmm = g{coeff_num};
        num_gaussians = numel(gmm.mu);

        % --- get present Gaussians from polyID_M ---
        if coeff_num <= numel(polyID_M) && ~isempty(polyID_M{coeff_num})
            gauss_indices = polyID_M{coeff_num};
        else
            % No relevant Gaussians → skip this coefficient
            continue;
        end

        % Clip invalid indices (must be within 1:num_gaussians)
        gauss_indices = gauss_indices(gauss_indices > 0 & gauss_indices <= num_gaussians);
        if isempty(gauss_indices)
            continue;  % skip coefficient if no valid Gaussians
        end

        % --- Loop over present Gaussians ---
        cluster_gauss_pct = zeros(num_clusters, num_gaussians);  % temp matrix
        for ii = 1:numel(gauss_indices)
            g_idx = gauss_indices(ii);  % scalar
            mu_g = gmm.mu(g_idx);
            S = gmm.Sigma;

            if ndims(S) == 3
                sigma_g = sqrt(squeeze(S(:,:,g_idx)));
            else
                sigma_g = sqrt(squeeze(S(g_idx,:)));
            end

            spike_vals = wd_coeff(:, coeff_num);
            range_st = mu_g - 2*abs(sigma_g);
            range_end = mu_g + 2*abs(sigma_g);
            in_range = (spike_vals >= range_st & spike_vals <= range_end);

            spike_clusters = cluster_times(in_range,1);

            % compute per-cluster %
            for c = 1:num_clusters
                cluster_id = all_clusters(c);
                n_in = sum(spike_clusters == cluster_id);
                if total_spikes_per_cluster(c) > 0
                    cluster_gauss_pct(c, g_idx) = (n_in / total_spikes_per_cluster(c)) * 100;
                end
            end
        end

        % --- take max percentage per cluster ---
        for c = 1:num_clusters
            [max_pct, max_idx_rel] = max(cluster_gauss_pct(c, gauss_indices), [], 2, 'omitnan');
            if max_pct > 0
                percentage_matrix(c, coeff_idx) = max_pct;
                winner_matrix(c, coeff_idx) = gauss_indices(max_idx_rel);
            else
                percentage_matrix(c, coeff_idx) = 0;
                winner_matrix(c, coeff_idx) = NaN;
            end
        end
    end

    % --- FIGURE 1: Percentage Heatmap ---
    fig1 = figure('Name', 'Cluster–Coefficient Max Percentage', 'NumberTitle', 'off', 'Position', [100 100 1200 700]);
    set(fig1, 'WindowStyle', 'docked');
    ax1 = axes('Parent', fig1);
    imagesc(ax1, percentage_matrix, [0 100]);
    colormap(ax1, turbo);
    colorbar(ax1);
    caxis(ax1, [0 100]);
    axis(ax1, 'tight', 'equal');

    xlabel(ax1, 'Coefficient');
    ylabel(ax1, 'Cluster');
    xticks(ax1, 1:num_coeffs);
    yticks(ax1, 1:num_clusters);
    xticklabels(ax1, arrayfun(@(x) sprintf('C%d', x), coeff_vals, 'UniformOutput', false));
    yticklabels(ax1, arrayfun(@(x) sprintf('Clust %d', x), all_clusters, 'UniformOutput', false));
    title(ax1, sprintf('Max Spike Percentage per Cluster–Coeff (Ch %s)', num2str(channelNum)), 'FontWeight','bold');

    % annotate values
    for c = 1:num_clusters
        for k = 1:num_coeffs
            val = percentage_matrix(c,k);
            if val > 0
                text_color = 'w';
                if val < 50, text_color = 'k'; end
                text(ax1, k, c, sprintf('%.1f', val), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'Color',text_color);
            end
        end
    end

    % --- FIGURE 2: Winner Gaussian Index Map ---
    fig2 = figure('Name', 'Cluster–Coefficient Gaussian Winner', 'NumberTitle', 'off', 'Position', [100 100 1200 700]);
    set(fig2, 'WindowStyle', 'docked');
    ax2 = axes('Parent', fig2);

    % compute max Gaussian index only from non-empty cells
    nonEmptyCells = ~cellfun(@isempty, polyID_M);
    max_gauss_idx = max(cellfun(@(x) max(x,[], 'omitnan'), polyID_M(nonEmptyCells)));

    cmap2 = lines(max_gauss_idx);  % distinct colors
    imagesc(ax2, winner_matrix, [1 max_gauss_idx]);
    colormap(ax2, cmap2);
    colorbar(ax2);
    axis(ax2, 'tight', 'equal');

    xlabel(ax2, 'Coefficient');
    ylabel(ax2, 'Cluster');
    xticks(ax2, 1:num_coeffs);
    yticks(ax2, 1:num_clusters);
    xticklabels(ax2, arrayfun(@(x) sprintf('C%d', x), coeff_vals, 'UniformOutput', false));
    yticklabels(ax2, arrayfun(@(x) sprintf('Clust %d', x), all_clusters, 'UniformOutput', false));
    title(ax2, sprintf('Winning Gaussian Index per Cluster–Coeff (Ch %s)', num2str(channelNum)), 'FontWeight','bold');

    % annotate Gaussian index numbers
    for c = 1:num_clusters
        for k = 1:num_coeffs
            g_idx = winner_matrix(c,k);
            if ~isnan(g_idx)
                text(ax2, k, c, sprintf('%d', g_idx), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'FontWeight','bold','Color','w');
            end
        end
    end

    % --- Save Outputs ---
    if ~exist(folderSpike, 'dir')
        mkdir(folderSpike);
    end
    % exportgraphics(fig1, fullfile(folderSpike, sprintf('ch%s_cluster_coeff_maxPct.png', num2str(channelNum))), 'Resolution', 300);
    % exportgraphics(fig2, fullfile(folderSpike, sprintf('ch%s_cluster_coeff_winner.png', num2str(channelNum))), 'Resolution', 300);

    fprintf('Saved two heatmaps to %s\n', folderSpike);
end
