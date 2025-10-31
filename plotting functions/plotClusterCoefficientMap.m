function plotClusterCoefficientMap(g, medDist_vec, wd_coeff, spikes, cluster_times, coeff_vals, select_gauss, ...
    channelNum,folderSpike,dataSetFlag)
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

    % Determine title/file suffix based on dataset
    switch dataSetFlag
        case 1
            dataStr = 'All Peaks';
            fileSuffix = 'allPks';
        case 2
            dataStr = 'Peaks Excluded';
            fileSuffix = 'excPks';
        otherwise
            error('dataSetFlag must be 1 (all peaks) or 2 (peaks excluded)');
    end

    %% format vecs to be cells
    % medDist_vec (or kDist_vec) assumed to have columns: [value, coeff#, gauss#]

    % Get unique coefficients in the data
    coeff_vals_data = unique(medDist_vec(:,2));
    num_coeffs_data = numel(coeff_vals_data);
    
    % Initialize polyID_M as cell array
    polyID_M = cell(max(coeff_vals_data),1);  % max coefficient number for indexing
    
    % Loop over each coefficient and collect Gaussian numbers
    for idx = 1:num_coeffs_data
        coeff = coeff_vals_data(idx);
        
        % Get unique Gaussians for this coefficient
        gauss_nums = unique(medDist_vec(medDist_vec(:,2) == coeff, 3));
        
        % Store in cell array at index = coeff number
        polyID_M{coeff} = gauss_nums;
    end
    
    % Now polyID_M_dyn can replace polyID_M in your function

    %%

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
            fprintf('Gauss empty: %d',coeff_num);
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

   % --- FIGURE 1: Percentage Heatmap (split by 16 coefficients per figure)
    maxCoeffsPerFig = 16;
    numFigures = ceil(num_coeffs / maxCoeffsPerFig);
    
    for fIdx = 1:numFigures
        coeffStart = (fIdx-1)*maxCoeffsPerFig + 1;
        coeffEnd   = min(fIdx*maxCoeffsPerFig, num_coeffs);
        coeffRange = coeffStart:coeffEnd;
        
        figF = figure('Name', sprintf('Cluster–Coefficient Max Percentage (Block %d)', fIdx), ...
                      'NumberTitle','off','Position',[100 100 1200 700]);
        set(figF, 'WindowStyle', 'docked');
        axF = axes('Parent', figF);
        
        % Extract sub-matrix
        pct_sub = percentage_matrix(:, coeffRange);
        
        imagesc(axF, pct_sub, [0 100]);
        colormap(axF, turbo);
        colorbar(axF);
        caxis(axF, [0 100]);
        axis(axF, 'tight', 'equal');
        
        xlabel(axF, 'Coefficient');
        ylabel(axF, 'Cluster');
        xticks(axF, 1:numel(coeffRange));
        yticks(axF, 1:num_clusters);
        xticklabels(axF, arrayfun(@(x) sprintf('C%d', coeff_vals(x)), coeffRange, 'UniformOutput', false));
        yticklabels(axF, arrayfun(@(x) sprintf('Clust %d', x), all_clusters, 'UniformOutput', false));
        title(axF, sprintf('Max Spike %% per Cluster–Coeff (%s, Ch %s)', dataStr, num2str(channelNum)), 'FontWeight','bold');
        
        % Annotate values
        for c = 1:num_clusters
            for k = 1:numel(coeffRange)
                val = pct_sub(c,k);
                if val > 0
                    text_color = 'w';
                    if val < 50, text_color = 'k'; end
                    text(axF, k, c, sprintf('%.1f', val), ...
                         'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'Color',text_color);
                end
            end
        end
        
        % Optional: save each block figure
        exportgraphics(figF, fullfile(folderSpike, sprintf('ch%s_cluster_coeff_maxPct_block%d_%s.png', num2str(channelNum), fIdx, fileSuffix)), 'Resolution',300);
    end
    % --- FIGURE 2: Winner Gaussian Index Map (split by 16 coefficients) ---
    maxCoeffsPerFig = 16;
    numFigures = ceil(num_coeffs / maxCoeffsPerFig);
    
    for fIdx = 1:numFigures
        coeffStart = (fIdx-1)*maxCoeffsPerFig + 1;
        coeffEnd   = min(fIdx*maxCoeffsPerFig, num_coeffs);
        coeffRange = coeffStart:coeffEnd;
    
        fig2 = figure('Name', sprintf('Cluster–Coefficient Gaussian Winner (Block %d)', fIdx), ...
                      'NumberTitle','off', 'Position', [100 100 1200 700]);
        set(fig2, 'WindowStyle', 'docked');
        ax2 = axes('Parent', fig2);
    
        % compute max Gaussian index only from non-empty cells
        nonEmptyCells = ~cellfun(@isempty, polyID_M);
        max_gauss_idx = max(cellfun(@(x) max(x,[], 'omitnan'), polyID_M(nonEmptyCells)));
    
        cmap2 = lines(max_gauss_idx);  % distinct colors
        imagesc(ax2, winner_matrix(:, coeffRange), [1 max_gauss_idx]);
        colormap(ax2, cmap2);
        colorbar(ax2);
        axis(ax2, 'tight', 'equal');
    
        xlabel(ax2, 'Coefficient');
        ylabel(ax2, 'Cluster');
        xticks(ax2, 1:numel(coeffRange));
        yticks(ax2, 1:num_clusters);
        xticklabels(ax2, arrayfun(@(x) sprintf('C%d', coeff_vals(x)), coeffRange, 'UniformOutput', false));
        yticklabels(ax2, arrayfun(@(x) sprintf('Clust %d', x), all_clusters, 'UniformOutput', false));
        title(ax2, sprintf('Winning Gaussian Index per Cluster–Coeff (%s, Ch %s)', dataStr, num2str(channelNum)), 'FontWeight','bold');
    
        % annotate Gaussian index numbers + spike counts
        for c = 1:num_clusters
            for k = 1:numel(coeffRange)
                coeff_num = coeff_vals(coeffRange(k));
                g_idx = winner_matrix(c, coeffRange(k));
                if ~isnan(g_idx)
                    % Plot Gaussian number
                    text(ax2, k, c, sprintf('%d', g_idx), ...
                        'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',8, 'FontWeight','bold', 'Color','w');
    
                    % Compute spikes in this gaussian for this cluster
                    spike_vals = wd_coeff(:, coeff_num);
                    gmm = g{coeff_num};
                    mu_g = gmm.mu(g_idx);
                    S = gmm.Sigma;
                    if ndims(S) == 3
                        sigma_g = sqrt(squeeze(S(:,:,g_idx)));
                    else
                        sigma_g = sqrt(squeeze(S(g_idx,:)));
                    end
                    sigma_g = abs(sigma_g(1));
                    in_gauss = spike_vals >= (mu_g - 2*sigma_g) & spike_vals <= (mu_g + 2*sigma_g);
                    cluster_mask = cluster_times(:,1) == all_clusters(c);
                    spikes_in_cluster = find(in_gauss & cluster_mask);
                    spikes_total = find(in_gauss);
    
                    % Add text below Gaussian number: selected / total spikes
                    txt_info = sprintf('%d/%d', numel(spikes_in_cluster), numel(spikes_total));
                    text(ax2, k, c+0.25, txt_info, ...
                        'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize',7, 'Color','w');
                end
            end
        end
    
        % --- Save figure ---
        if ~exist(folderSpike, 'dir')
            mkdir(folderSpike);
        end
        exportgraphics(fig2, fullfile(folderSpike, sprintf('ch%s_cluster_coeff_winner_block%d_%s.png', num2str(channelNum), fIdx, fileSuffix)), 'Resolution',300);
        fprintf('Saved figure block %d to %s\n', fIdx, folderSpike);
    end

    % select_gauss assumed to be [medDist, k_value, coeff#, gauss#]
    
    % Extract unique coefficients and gaussians from select_gauss
    selected_coeffs = unique(select_gauss(:,3));
    selected_gauss  = select_gauss(:,4);
    
    % Initialize filtered matrices
    percentage_matrix_sel = zeros(num_clusters, numel(selected_coeffs));
    winner_matrix_sel     = nan(num_clusters, numel(selected_coeffs));
    

       for idx = 1:numel(selected_coeffs)
            coeff = selected_coeffs(idx);
        
            % find gaussian rows for this coeff in select_gauss and get unique gaussians
            rows = find(select_gauss(:,3) == coeff);
            if isempty(rows)
                continue;
            end
            gauss_list = unique(select_gauss(rows,4));
        
            % skip if no GMM for this coeff
            if coeff > numel(g) || isempty(g{coeff}) || isempty(g{coeff}.mu)
                continue;
            end
            gmm = g{coeff};
            num_gauss_total = numel(gmm.mu);
        
            % clip gauss_list to valid indices
            gauss_list = gauss_list(gauss_list>0 & gauss_list<=num_gauss_total);
            if isempty(gauss_list)
                continue;
            end
        
            % preallocate per-gaussian per-cluster percent matrix (num_clusters x num_selected_gauss)
            cluster_gauss_pct = zeros(num_clusters, numel(gauss_list));
        
            spike_vals = wd_coeff(:, coeff);  % values for this coefficient (all spikes)
        
            for gg = 1:numel(gauss_list)
                g_idx = gauss_list(gg);
                mu_g = gmm.mu(g_idx);
                S = gmm.Sigma;
        
                % compute sigma robustly (handles diag/3D)
                if ndims(S) == 3
                    sigma_g = sqrt(squeeze(S(:,:,g_idx)));
                else
                    % S might be stored per component as a row
                    sigma_g = sqrt(squeeze(S(g_idx,:)));
                end
                sigma_g = abs(sigma_g);  % ensure positive
        
                % define selection range (2*sigma either scalar or vector; use scalar magnitude)
                % If sigma_g is vector, take scalar std (mean or first element). Use first element.
                if numel(sigma_g) > 1
                    sg = sigma_g(1);
                else
                    sg = sigma_g;
                end
                range_st = mu_g - 2*sg;
                range_end = mu_g + 2*sg;
        
                in_range = (spike_vals >= range_st & spike_vals <= range_end);
                if ~any(in_range)
                    continue;
                end
        
                spike_clusters = cluster_times(in_range,1);  % cluster IDs of spikes in this gaussian
        
                % compute counts per cluster in the same order as all_clusters
                for c = 1:num_clusters
                    cluster_id = all_clusters(c);
                    cluster_gauss_pct(c, gg) = 0;  % default zero
                    if total_spikes_per_cluster(c) == 0
                        continue;
                    end
                    n_in = sum(spike_clusters == cluster_id);
                    cluster_gauss_pct(c, gg) = (n_in / total_spikes_per_cluster(c)) * 100;
                end
            end
        
            % Now pick max percentage across selected gaussians for each cluster
            % cluster_gauss_pct is num_clusters x num_selected_gauss
            [max_pct_vec, max_idx_rel] = max(cluster_gauss_pct, [], 2, 'omitnan');  % column index of best gaussian in gauss_list
        
            % populate selection matrices
            for c = 1:num_clusters
                if max_pct_vec(c) > 0
                    percentage_matrix_sel(c, idx) = max_pct_vec(c);
                    winner_matrix_sel(c, idx) = gauss_list(max_idx_rel(c));
                else
                    percentage_matrix_sel(c, idx) = 0;
                    winner_matrix_sel(c, idx) = NaN;
                end
            end
        end

    
    % --- FIGURE 3: Percentage Heatmap for Selected Gaussians (split by 16 coefficients) ---
    maxCoeffsPerFig = 16;
    numSelCoeffs = numel(selected_coeffs);
    numFigures = ceil(numSelCoeffs / maxCoeffsPerFig);
    
    for fIdx = 1:numFigures
        coeffStart = (fIdx-1)*maxCoeffsPerFig + 1;
        coeffEnd   = min(fIdx*maxCoeffsPerFig, numSelCoeffs);
        coeffRange = coeffStart:coeffEnd;
    
        figSel = figure('Name', sprintf('Selected Gaussians Percentage (Block %d)', fIdx), ...
                        'NumberTitle','off', 'Position',[100 100 1200 700]);
        set(figSel, 'WindowStyle', 'docked');
        axSel = axes('Parent', figSel);
    
        % extract submatrix
        pct_sub = percentage_matrix_sel(:, coeffRange);
    
        imagesc(axSel, pct_sub, [0 100]);
        colormap(axSel, turbo);
        colorbar(axSel);
        caxis(axSel, [0 100]);
        axis(axSel, 'tight', 'equal');
    
        xlabel(axSel,'Selected Coefficients');
        ylabel(axSel,'Cluster');
        xticks(axSel, 1:numel(coeffRange));
        yticks(axSel, 1:num_clusters);
        xticklabels(axSel, arrayfun(@(x) sprintf('C%d', selected_coeffs(x)), coeffRange, 'UniformOutput', false));
        yticklabels(axSel, arrayfun(@(x) sprintf('Clust %d', x), all_clusters, 'UniformOutput', false));
        title(axSel, sprintf('Max Spike %% for Selected Gaussians (%s, Ch %s)', dataStr, num2str(channelNum)), 'FontWeight','bold');
    
        % annotate values
        for c = 1:num_clusters
            for k = 1:numel(coeffRange)
                val = pct_sub(c,k);
                if val > 0
                    text(axSel, k, c, sprintf('%.1f', val), ...
                         'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'Color','k');
                end
            end
        end
    
        % save each block
        if ~exist(folderSpike, 'dir')
            mkdir(folderSpike);
        end
        exportgraphics(figSel, fullfile(folderSpike, sprintf('ch%s_selectedGauss_pct_block%d_%s.png', num2str(channelNum), fIdx, fileSuffix)), 'Resolution',300);
        fprintf('Saved selected Gauss percentage figure block %d to %s\n', fIdx, folderSpike);
    end

    % --- FIGURE 4: Winner Gaussian Index for selected Gaussians (segmented + annotated) ---
    maxCoeffsPerFig = 16;
    numSelCoeffs = numel(selected_coeffs);
    numFigures = ceil(numSelCoeffs / maxCoeffsPerFig);
    
    for fIdx = 1:numFigures
        coeffStart = (fIdx-1)*maxCoeffsPerFig + 1;
        coeffEnd   = min(fIdx*maxCoeffsPerFig, numSelCoeffs);
        coeffRange = coeffStart:coeffEnd;
    
        figWinnerSel = figure('Name', sprintf('Selected Gaussians Winner (Block %d)', fIdx), ...
                              'NumberTitle','off','Position',[100 100 1200 700]);
        set(figWinnerSel, 'WindowStyle', 'docked');
        axWinnerSel = axes('Parent', figWinnerSel);
    
        % extract submatrix
        winner_sub = winner_matrix_sel(:, coeffRange);
    
        % colormap
        max_gauss_sel = max(select_gauss(:,4));
        cmapSel = lines(max_gauss_sel);
        imagesc(axWinnerSel, winner_sub, [1 max_gauss_sel]);
        colormap(axWinnerSel, cmapSel);
        colorbar(axWinnerSel);
        axis(axWinnerSel, 'tight', 'equal');
    
        xlabel(axWinnerSel,'Selected Coefficients');
        ylabel(axWinnerSel,'Cluster');
        xticks(axWinnerSel, 1:numel(coeffRange));
        yticks(axWinnerSel, 1:num_clusters);
        xticklabels(axWinnerSel, arrayfun(@(x) sprintf('C%d', selected_coeffs(x)), coeffRange, 'UniformOutput', false));
        yticklabels(axWinnerSel, arrayfun(@(x) sprintf('Clust %d', x), all_clusters, 'UniformOutput', false));
        title(axWinnerSel, sprintf('Winning Gaussian per Cluster (Block %d, %s, Ch %s)', fIdx, dataStr, num2str(channelNum)), 'FontWeight','bold');
    
        % --- Annotate each cell ---
        for c = 1:num_clusters
            for k = 1:numel(coeffRange)
                coeff_pos = coeffRange(k);          % position in selected_coeffs
                coeff = selected_coeffs(coeff_pos); % actual coefficient number
                g_idx = winner_matrix_sel(c, coeff_pos);
                if ~isnan(g_idx)
                    gmm = g{coeff};
                    S = gmm.Sigma;
                    if ndims(S) == 3
                        sigma_g = sqrt(squeeze(S(:,:,g_idx)));
                    else
                        sigma_g = sqrt(squeeze(S(g_idx,:)));
                    end
    
                    mu_g = gmm.mu(g_idx);
                    range_st = mu_g - 2*sigma_g;
                    range_end = mu_g + 2*sigma_g;
    
                    spike_vals = wd_coeff(:, coeff);
    
                    spikes_in_cluster = sum((spike_vals >= range_st & spike_vals <= range_end) & cluster_times(:,1) == all_clusters(c));
                    total_spikes_gauss = sum(spike_vals >= range_st & spike_vals <= range_end);
    
                    text_str = sprintf('%d\n[%d/%d]', g_idx, spikes_in_cluster, total_spikes_gauss);
                    text(axWinnerSel, k, c, text_str, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',8, 'FontWeight','bold','Color','w');
                end
            end
        end
    
        % --- Save each block ---
        if ~exist(folderSpike, 'dir')
            mkdir(folderSpike);
        end
        exportgraphics(figWinnerSel, fullfile(folderSpike, sprintf('ch%s_selectedGauss_winner_block%d_%s.png', num2str(channelNum), fIdx, fileSuffix)), 'Resolution',300);
        fprintf('Saved winner heatmap block %d to %s\n', fIdx, folderSpike);
    end


end
