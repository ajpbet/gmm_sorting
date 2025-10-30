function tab_gauss(summary_table, channelNum, coeff_idx, g, folderName)    % TAB_GAUSS: Visualize Gaussian coefficients as color-coded tables (global + local scale)
    % with weight filtering and M_comp separation
    %
    % Inputs:
    %   summary_table.kj   - Cell array where each cell has 8 numeric values
    %   channelNum (opt)   - Channel number for title/label (default = 1)
    %   coeff_idx (opt)    - Vector of coefficient indices to plot (default = all)
    %   g (opt)            - GMM objects for weight filtering

    if nargin < 2 || isempty(channelNum)
        channelNum = 1;
    end

    if nargin < 3 || isempty(coeff_idx)
        coeff_idx = 1:numel(summary_table.kj);
    end

    if nargin < 4 || isempty(g)
        % If no GMM objects provided, use all components
        weight_filter = false;
    else
        weight_filter = true;
    end

    if nargin < 5 || isempty(folderName)
        folderName = pwd;  % Default to current directory
    end

    % Create directory if it doesn't exist
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end

    kj_all = summary_table.kj(coeff_idx);
    numCoeff = numel(kj_all);

    % Apply weight filtering and separate M_comp vs non-M_comp
    [gaussVals_mcomp, gaussVals_nonmcomp, globalMin, globalMax] = filter_and_separate_data(summary_table, kj_all, g, weight_filter, coeff_idx);

    baseWidth = 300;
    figWidth = max(800, baseWidth + 25 * numCoeff);

    % -------- Global colormap - M_comp --------
    % if ~isempty(gaussVals_mcomp)
    %     fig1 = figure('Name', sprintf('Channel %s - M_comp Gaussian Coefficients (Global Color Scale)', channelNum), ...
    %                   'Color', 'w', 'Position', [100, 300, figWidth, 500], 'WindowStyle', 'docked');
    %     ax1 = axes('Parent', fig1);
    %     imagesc(ax1, gaussVals_mcomp);
    %     colormap(ax1, "parula");
    %     colorbar;
    %     caxis([globalMin, globalMax]);
    % 
    %     title(ax1, sprintf('Channel %s: M_comp Gaussian Parameters (Global Scale) [%d Coefficients, %d Components]', ...
    %           channelNum, numCoeff, size(gaussVals_mcomp, 1)), 'FontSize', 12, 'FontWeight', 'bold');
    %     format_axes(ax1, summary_table, coeff_idx, numCoeff, size(gaussVals_mcomp, 1), 12);
    %     add_labels(ax1, gaussVals_mcomp, 10);
    % 
    %     exportgraphics(fig1, fullfile(folderName,sprintf('Channel%s_Mcomp_Gaussian_Global.png', channelNum)), 'Resolution', 300);
    % end

    % -------- Global colormap - Non-M_comp --------
    if ~isempty(gaussVals_nonmcomp)
        fig2 = figure('Name', sprintf('Channel %s - Non-M_comp Gaussian Coefficients (Global Color Scale)', channelNum), ...
                      'Color', 'w', 'Position', [150, 300, figWidth, 500], 'WindowStyle', 'docked');
        ax2 = axes('Parent', fig2);
        imagesc(ax2, gaussVals_nonmcomp);
        colormap(ax2, "parula");
        colorbar;
        caxis([globalMin, globalMax]);

        title(ax2, sprintf('Channel %s: Non-M_comp Gaussian Parameters (Global Scale) [%d Coefficients, %d Components]', ...
              channelNum, numCoeff, size(gaussVals_nonmcomp, 1)), 'FontSize', 12, 'FontWeight', 'bold');
        format_axes(ax2, summary_table, coeff_idx, numCoeff, size(gaussVals_nonmcomp, 1), 12);
        add_labels(ax2, gaussVals_nonmcomp, 10);

        exportgraphics(fig2, fullfile(folderName,sprintf('Channel%s_NonMcomp_Gaussian_Global.png', channelNum)), 'Resolution', 300);
    end

    % -------- Local colormap - M_comp --------
    % if ~isempty(gaussVals_mcomp)
    %     localVals_mcomp = normalize_local(gaussVals_mcomp);
    % 
    %     fig3 = figure('Name', sprintf('Channel %s - M_comp Gaussian Coefficients (Local Color Scale)', channelNum), ...
    %                   'Color', 'w', 'Position', [200, 300, figWidth, 500], 'WindowStyle', 'docked');
    %     ax3 = axes('Parent', fig3);
    %     imagesc(ax3, localVals_mcomp);
    %     colormap(ax3,"parula");
    %     colorbar;
    %     caxis([0, 1]);
    % 
    %     title(ax3, sprintf('Channel %s: M_comp Gaussian Parameters (Local Scale) [%d Coefficients, %d Components]', ...
    %           channelNum, numCoeff, size(gaussVals_mcomp, 1)), 'FontSize', 12, 'FontWeight', 'bold');
    %     format_axes(ax3, summary_table, coeff_idx, numCoeff, size(gaussVals_mcomp, 1), 12);
    %     add_labels(ax3, gaussVals_mcomp, 10);
    % 
    %     exportgraphics(fig3, fullfile(folderName,sprintf('Channel%s_Mcomp_Gaussian_Local.png', channelNum)), 'Resolution', 300);
    % end

    % -------- Local colormap - Non-M_comp --------
    if ~isempty(gaussVals_nonmcomp)
        localVals_nonmcomp = normalize_local(gaussVals_nonmcomp);
        
        fig4 = figure('Name', sprintf('Channel %s - Non-M_comp Gaussian Coefficients (Local Color Scale)', channelNum), ...
                      'Color', 'w', 'Position', [250, 300, figWidth, 500], 'WindowStyle', 'docked');
        ax4 = axes('Parent', fig4);
        imagesc(ax4, localVals_nonmcomp);
        colormap(ax4, "parula");
        colorbar;
        caxis([0, 1]);

        title(ax4, sprintf('Channel %s: Non-M_comp Gaussian Parameters (Local Scale) [%d Coefficients, %d Components]', ...
              channelNum, numCoeff, size(gaussVals_nonmcomp, 1)), 'FontSize', 12, 'FontWeight', 'bold');
        format_axes(ax4, summary_table, coeff_idx, numCoeff, size(gaussVals_nonmcomp, 1), 12);
        add_labels(ax4, gaussVals_nonmcomp, 10);

        exportgraphics(fig4, fullfile(folderName,sprintf('Channel%s_NonMcomp_Gaussian_Local.png', channelNum)), 'Resolution', 300);
    end
end

% --- Helper function for weight filtering and M_comp separation ---
function [gaussVals_mcomp, gaussVals_nonmcomp, globalMin, globalMax] = filter_and_separate_data(summary_table, kj_all, g, weight_filter, coeff_idx)
    numCoeff = numel(kj_all);
    all_filtered_values = [];
    
    % Pre-calculate which components to keep for each coefficient
    filtered_data = cell(numCoeff, 1);
    filtered_positions_all = cell(numCoeff, 1);
    
    for i = 1:numCoeff
        current_kj = kj_all{i};
        current_M_comp = summary_table.M_comp{i};
        
        if weight_filter && ~isempty(g{i})
            % Apply weight filtering
            g_in = g{i};
            comp_proportions = g_in.ComponentProportion;
            keep = comp_proportions > 0.005;
            filtered_kj = current_kj(keep);
            filtered_positions = 1:length(current_kj);
            filtered_positions = filtered_positions(keep);
        else
            % No filtering
            filtered_kj = current_kj;
            filtered_positions = 1:length(current_kj);
        end
        
        filtered_data{i} = filtered_kj;
        filtered_positions_all{i} = filtered_positions;
        all_filtered_values = [all_filtered_values; filtered_kj(:)];
    end
    
    globalMin = min(all_filtered_values);
    globalMax = max(all_filtered_values);
    
    % Find all unique component indices across all coefficients
    all_component_indices = unique(cell2mat(cellfun(@(x) x(:), filtered_positions_all, 'UniformOutput', false)));
    
    % Initialize arrays with NaN (so we can handle missing values)
    gaussVals_mcomp = nan(length(all_component_indices), numCoeff);
    gaussVals_nonmcomp = nan(length(all_component_indices), numCoeff);
    
    for i = 1:numCoeff
        current_kj = filtered_data{i};
        current_positions = filtered_positions_all{i};
        current_M_comp = summary_table.M_comp{i};
        
        % Find where these positions exist in our master component list
        [~, idx_in_master] = ismember(current_positions, all_component_indices);
        
        % Separate M_comp vs non-M_comp
        is_mcomp = ismember(current_positions, current_M_comp);
        
        % Assign values to the appropriate arrays
        for j = 1:length(current_positions)
            master_idx = idx_in_master(j);
            if is_mcomp(j)
                gaussVals_mcomp(master_idx, i) = current_kj(j);
            else
                gaussVals_nonmcomp(master_idx, i) = current_kj(j);
            end
        end
    end
    
    % Remove rows that are all NaN (components that don't appear in any coefficient)
    gaussVals_mcomp = gaussVals_mcomp(~all(isnan(gaussVals_mcomp), 2), :);
    gaussVals_nonmcomp = gaussVals_nonmcomp(~all(isnan(gaussVals_nonmcomp), 2), :);
end
% --- Helper function for local normalization ---
function localVals = normalize_local(values)
    localVals = zeros(size(values));
    for j = 1:size(values, 2)
        col = values(:, j);
        if range(col) > 0
            localVals(:, j) = (col - min(col)) / (max(col) - min(col));
        else
            localVals(:, j) = 0.5;
        end
    end
end

% --- Existing helper functions (keep the same) ---
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

% function cmap = redgreencmap()
%     n = 256;
%     % More intense/vibrant version
%     mid_point = round(n/2);
% 
%     % Deep blue to cyan to white
%     r_neg = linspace(0, 0.2, mid_point)';
%     g_neg = linspace(0, 0.8, mid_point)';
%     b_neg = linspace(0.8, 1, mid_point)';
% 
%     % White to yellow to deep red
%     r_pos = linspace(1, 0.9, n-mid_point)';
%     g_pos = linspace(1, 0.1, n-mid_point)';
%     b_pos = linspace(1, 0.1, n-mid_point)';
% 
%     cmap = [r_neg, g_neg, b_neg; r_pos, g_pos, b_pos];
% end