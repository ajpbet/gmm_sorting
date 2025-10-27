function tab_gauss(summary_table, channelNum, coeff_idx)
    % TAB_GAUSS: Visualize Gaussian coefficients as color-coded tables (global + local scale)
    %
    % Inputs:
    %   summary_table.kj   - Cell array where each cell has 8 numeric values
    %   channelNum (opt)   - Channel number for title/label (default = 1)
    %   coeff_idx (opt)    - Vector of coefficient indices to plot (default = all)

    if nargin < 2 || isempty(channelNum)
        channelNum = 1;
    end

    if nargin < 3 || isempty(coeff_idx)
        coeff_idx = 1:numel(summary_table.kj);
    end

    kj_all = summary_table.kj(coeff_idx);
    numCoeff = numel(kj_all);

    try
        gaussVals = cell2mat(cellfun(@(x) x(:), kj_all(:)', 'UniformOutput', false));
    catch
        error('Each summary_table.kj cell must contain exactly 8 numeric values.');
    end
    numGauss = size(gaussVals, 1);

    globalMin = min(cell2mat(summary_table.kj(:)));
    globalMax = max(cell2mat(summary_table.kj(:)));

    baseWidth = 300;
    figWidth = max(800, baseWidth + 25 * numCoeff);

    % -------- Global colormap --------
    fig1 = figure('Name', sprintf('Channel %d - Gaussian Coefficients (Global Color Scale)', channelNum), ...
                  'Color', 'w', 'Position', [100, 300, figWidth, 500], 'WindowStyle', 'docked');
    ax1 = axes('Parent', fig1);
    imagesc(ax1, gaussVals);
    colormap(ax1, redgreencmap());
    colorbar;
    caxis([globalMin, globalMax]);

    title(ax1, sprintf('Channel %d: Gaussian Parameters (Global Scale) [%d Coefficients]', channelNum, numCoeff), ...
        'FontSize', 12, 'FontWeight', 'bold');
    format_axes(ax1, summary_table, coeff_idx, numCoeff, numGauss, 12);
    add_labels(ax1, gaussVals, 10);

    % Save docked figure with correct size
    exportgraphics(fig1, sprintf('Channel%d_Gaussian_Global.png', channelNum), 'Resolution', 300);

    % -------- Local colormap --------
    localVals = zeros(size(gaussVals));
    for j = 1:numCoeff
        col = gaussVals(:, j);
        if range(col) > 0
            localVals(:, j) = (col - min(col)) / (max(col) - min(col));
        else
            localVals(:, j) = 0.5;
        end
    end

    fig2 = figure('Name', sprintf('Channel %d - Gaussian Coefficients (Local Color Scale)', channelNum), ...
                  'Color', 'w', 'Position', [150, 300, figWidth, 500], 'WindowStyle', 'docked');
    ax2 = axes('Parent', fig2);
    imagesc(ax2, localVals);
    colormap(ax2, redgreencmap());
    colorbar;
    caxis([0, 1]);

    title(ax2, sprintf('Channel %d: Gaussian Parameters (Local Scale) [%d Coefficients]', channelNum, numCoeff), ...
        'FontSize', 12, 'FontWeight', 'bold');
    format_axes(ax2, summary_table, coeff_idx, numCoeff, numGauss, 12);
    add_labels(ax2, gaussVals, 10);

    % Save docked figure with correct size
    exportgraphics(fig2, sprintf('Channel%d_Gaussian_Local.png', channelNum), 'Resolution', 300);
end

% --- helper functions ---
function format_axes(ax, summary_table, coeff_idx, numCoeff, numGauss, fontSize)
    xlabel(ax, 'Coefficients', 'FontSize', fontSize);
    ylabel(ax, 'Gaussians', 'FontSize', fontSize);
    xticks(ax, 1:numCoeff);
    if isfield(summary_table, 'coeff_labels')
        xticklabels(ax, summary_table.coeff_labels(coeff_idx));
    else
        xticklabels(ax, arrayfun(@(i) sprintf('Coeff%d', coeff_idx(i)), 1:numCoeff, 'UniformOutput', false));
    end
    yticks(ax, 1:numGauss);
    yticklabels(ax, arrayfun(@(i) sprintf('Gauss%d', i), 1:numGauss, 'UniformOutput', false));
    if numCoeff > 10
        xtickangle(ax, 45);
    end
    set(ax, 'FontSize', fontSize, 'Box', 'on');
end

function add_labels(ax, values, fontSize)
    [numGauss, numCoeff] = size(values);
    hold(ax, 'on');
    for i = 1:numGauss
        for j = 1:numCoeff
            val = values(i, j);
            text(ax, j, i, sprintf('%.2f', val), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', fontSize, ...
                'Color', 'k', ...
                'FontWeight', 'bold');
        end
    end
    hold(ax, 'off');
end

function cmap = redgreencmap()
    n = 256;
    r = linspace(1, 0, n)';
    g = linspace(0, 1, n)';
    b = zeros(n, 1);
    cmap = [r g b];
end
