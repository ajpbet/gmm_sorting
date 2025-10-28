function plotGaussianClusterAnalysis(g, polyID_M, wd_coeff, cluster_times, coeff_vals, channelNum, folderSpike,M_comp)
% PLOTGAUSSIANCLUSTERANALYSIS
% Visualizes Gaussian vs Coefficient dominance.
% Rows = original Gaussian IDs (from polyID_M)
% Columns = coefficients
% Cell color = dominant cluster color, brightness = normalized percent.

    % Cluster colors (cluster 0 = first row)
    clus_colors = [ ...
        0.0 0.0 0.0; ...
        0.0 0.0 1.0; ...
        1.0 0.0 0.0; ...
        0.0 0.5 0.0; ...
        0.6 0.0 0.0; ...
        0.4 0.0 0.75; ...
        0.96 0.52 0.03; ...
        0.45 0.38 0.24; ...
        1.0 0.1 0.72; ...
        0.55 0.55 0.55 ...
    ];

    all_clusters = unique(cluster_times(:,1));
    total_per_cluster = accumarray(cluster_times(:,1)+1, 1);
    num_coeffs = numel(coeff_vals);

    % Determine global Gaussian IDs
    all_gauss_ids = unique(cell2mat(polyID_M(:)));
    num_gauss_total = max(all_gauss_ids);

    bestClusMat = nan(num_gauss_total, num_coeffs);
    bestPctMat  = nan(num_gauss_total, num_coeffs);

    % === Compute dominant cluster for each Gaussian ===
    for c = 1:num_coeffs
        coeff_num = coeff_vals(c);
        if isempty(g{coeff_num}) || isempty(polyID_M{coeff_num}), continue; end

        mu = g{coeff_num}.mu(:);
        sigma = sqrt(squeeze(g{coeff_num}.Sigma(:)));
        gaussIDs = polyID_M{coeff_num};

        for k = 1:numel(mu)
            g_id = gaussIDs(k); % original Gaussian number
            range = [mu(k)-2*sigma(k), mu(k)+2*sigma(k)];
            idx = (wd_coeff(:,coeff_num) >= range(1)) & (wd_coeff(:,coeff_num) <= range(2));
            if ~any(idx)||ismember(g_id,M_comp{c}), continue; end

            clusIDs = cluster_times(idx,1);
            counts = accumarray(clusIDs+1, 1, [max(all_clusters)+1, 1]);
            perc = (counts ./ total_per_cluster) * 100;

            [max_pct, best_clus] = max(perc);
            if max_pct > 0
                bestClusMat(g_id, c) = best_clus - 1;
                bestPctMat(g_id, c)  = max_pct;
            end
        end
    end

    % === Normalize for brightness scaling ===
    global_max = max(bestPctMat(:), [], 'omitnan');
    if isempty(global_max) || global_max == 0, global_max = 1; end

    % === Plot table ===
    fig = figure('Name', 'Gaussian Cluster Table', ...
                 'NumberTitle', 'off', 'Color', 'w', ...
                 'Position', [100 100 1600 900]);
    ax = axes(fig);
    hold(ax, 'on');
    axis(ax, 'equal');

    [nG, nC] = size(bestClusMat);
    square_size = 0.9;

    for c = 1:nC
        for g_id = 1:nG
            clus = bestClusMat(g_id, c);
            pct  = bestPctMat(g_id, c);
            if isnan(clus) || isnan(pct), continue; end

            color_idx = min(clus+1, size(clus_colors,1));
            base_col = clus_colors(color_idx,:);
            bright = pct / global_max;
            color = base_col * bright + (1-bright)*[1 1 1];

            rectangle('Position',[c-0.5, g_id-0.5, square_size, square_size], ...
                      'FaceColor', color, 'EdgeColor', 'k');
            text_color = 'k';
            if bright > 0.6, text_color = 'w'; end
            text(c, g_id, sprintf('%d\n%.1f%%', clus, pct), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color',text_color);
        end
    end

    % === Axes and labels ===
    xlim([0.5, nC+0.5]);
    ylim([0.5, nG+0.5]);
    xticks(1:nC);
    xticklabels(arrayfun(@(x) sprintf('C%d', x), coeff_vals, 'UniformOutput', false));
    yticks(1:nG);
    yticklabels(arrayfun(@(x) sprintf('G%d', x), 1:nG, 'UniformOutput', false));
    xlabel('Coefficient Number');
    ylabel('Gaussian Number');
    title(sprintf('Gaussian Cluster Dominance (Channel %s)\nColor = Cluster, Brightness = Relative %%', channelNum), ...
          'FontWeight', 'bold');

    box on;
    hold off;

    % === Legend below the table (outside axes) ===
    num_clusters = min(size(clus_colors,1), max(all_clusters)+1);
    leg_box_w = 0.03;  % width of each color box (normalized units)
    leg_box_h = 0.03;  % height of each color box
    leg_gap_x = 0.01;  % horizontal gap between boxes

    % Make figure a bit taller to fit legend
    figPos = get(fig, 'Position');
    figPos(4) = figPos(4) * 1.2; % 20% taller
    set(fig, 'Position', figPos);

    % Move axes up slightly to make space
    axPos = get(ax, 'Position');
    axPos(2) = axPos(2) + 0.10;
    axPos(4) = axPos(4) * 0.85;
    set(ax, 'Position', axPos);

    % Compute total legend width (normalized)
    total_leg_width = num_clusters * (leg_box_w + leg_gap_x);
    start_x = 0.5 - total_leg_width / 2; % center it
    leg_y = 0.05; % near bottom of figure

    % Draw each cluster box and label
    for i = 1:num_clusters
        x_rel = start_x + (i-1) * (leg_box_w + leg_gap_x);
        annotation(fig, 'rectangle', [x_rel, leg_y, leg_box_w, leg_box_h], ...
                   'FaceColor', clus_colors(i,:), 'EdgeColor', 'k');
        annotation(fig, 'textbox', [x_rel, leg_y - 0.035, leg_box_w, 0.03], ...
                   'String', sprintf('%d', i-1), ...
                   'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                   'VerticalAlignment', 'top', 'FontSize', 9, 'FontWeight', 'bold');
    end

    hold off;

    % === Save both figures ===
    if ~exist(folderSpike, 'dir')
        mkdir(folderSpike);
    end
    mainFig = fullfile(folderSpike, sprintf('ch%s_gaussClusterTable.png', channelNum));
    legFig  = fullfile(folderSpike, sprintf('ch%s_clusterLegend.png', channelNum));
    exportgraphics(fig, mainFig, 'Resolution', 300);

    fprintf('Saved:\n- %s\n- %s\n', mainFig, legFig);
end
