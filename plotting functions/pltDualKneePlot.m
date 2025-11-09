function pltDualKneePlot(len1, vector1, len2, vector2, ks_coeff, select_gauss1, select_gauss2, plotID1, plotID2, channelNum, folderName)
% PLTDUALKNEEPLOT plots two vectors on a dual Y-axis figure, applying a
% selection logic based on KS coefficients and dual Gaussian selections.
%
% Inputs:
%   len1:           Length of the 'selected' portion of vector1 (for cutoff line).
%   vector1:        Full sorted vector 1: [value, coeff#, gauss#]
%   len2:           Length of the 'selected' portion of vector2 (for cutoff line).
%   vector2:        Full sorted vector 2: [value, coeff#, gauss#]
%   ks_coeff:       Vector of coefficient indices selected by KS test.
%   select_gauss1:  Matrix 1: [kv, meddist, coeff#, gauss#] - first selection boundary
%   select_gauss2:  Matrix 2: [kv, meddist, coeff#, gauss#] - second selection boundary (subset)
%   plotID1:        String identifier for vector 1 (e.g., 'kv').
%   plotID2:        String identifier for vector 2 (e.g., 'medDist').
%   channelNum:     Channel number for filename.
%   folderName:     Folder to save the output file.

    % Define original colors (for lines and default classification)
    COLOR_1 = [0.0 0.447 0.741]; % MATLAB Blue for Left Axis
    COLOR_2 = [0.85 0.325 0.098]; % MATLAB Red/Orange for Right Axis
    
    % Define new colors for classification categories
    COLOR_KS = [0 0 0];          % Black for KS coefficients (Triangle marker)
    COLOR_SG1_AND_SG2 = [0 0.6 0]; % Dark Green (Inside both SG1 and SG2 boundaries)
    COLOR_SG1_ONLY = [0.8 0 0.8];  % Magenta (Inside SG1 but NOT SG2 boundary)

    fig = figure('Visible', 'on');
    hold on;
    grid on;
    
    title(sprintf('%s vs %s per Coefficient (Channel %s)', plotID1, plotID2, num2str(channelNum)));
    xlabel('Coefficient Index (Sorted by Value)');

title(sprintf('%s vs %s per Coefficient (Channel %s)', plotID1, plotID2, num2str(channelNum)));
    xlabel('Coefficient Index (Sorted by Value)');
    
    % select_gauss matrices: Coeff indices in column 3, Gaussian indices in column 4
    count_coeff_sg1 = length(unique(select_gauss1(:, 3)));
    count_gauss_sg1 = length(select_gauss1(:, 4));
    count_coeff_sg2 = length(unique(select_gauss2(:, 3)));
    count_gauss_sg2 = length(select_gauss2(:, 4));
    count_ks = length(unique(ks_coeff));
    
    % Using 'figure' units to place text in the figure margin, outside the axes.
    text_x = 0.02; % Normalized X coordinate relative to FIGURE (right margin)
    text_y_start = 0.7; % Normalized Y coordinate relative to FIGURE
    line_spacing = 0.04;
    
    % SG1 Counts (Before removal)
    text(text_x, text_y_start, ...
        sprintf('Init Coeff Select: %d', count_coeff_sg1), ...
        'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'left');
    
    text(text_x, text_y_start - line_spacing, ...
        sprintf('Init Gauss Select: %d', count_gauss_sg1), ...
        'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'left');
        
    % SG2 Counts (After removal)
    text(text_x, text_y_start - 2*line_spacing, ...
        sprintf('Coeffs After Removal: %d', count_coeff_sg2), ...
        'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'left');
        
    text(text_x, text_y_start - 3*line_spacing, ...
        sprintf('Gauss After Removal: %d', count_gauss_sg2), ...
        'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'left');
        
    % KS Coeff count (Explicitly requested for the right side)
    text(text_x, text_y_start - 4*line_spacing, ...
        sprintf('KS Coeff Count: %d', count_ks), ...
        'Units', 'normalized', 'FontSize', 9, 'HorizontalAlignment', 'left');

    % Calculate cutoff points (used only for drawing vertical lines)
    cutoff1 = length(vector1) - len1;
    cutoff2 = length(vector2) - len2;

    
    % Get unique coefficient indices from the selection matrices (Col 3)
    SG1_Coeffs = select_gauss1(:, 3);
    SG2_Coeffs = select_gauss2(:, 3);
    
    V1_Coeffs = vector1(:, 2);
    is_ks_v1 = ismember(V1_Coeffs, ks_coeff);
    is_sg1_v1 = ismember(V1_Coeffs, SG1_Coeffs);
    is_sg2_v1 = ismember(V1_Coeffs, SG2_Coeffs);
    
    % Mutually exclusive masks for V1
    mask_ks_v1 = is_ks_v1;                                    % Group 1: KS (Highest priority)
    mask_sg1_only_v1 = is_sg1_v1 & ~is_sg2_v1 & ~is_ks_v1;     % Group 2: SG1 Only (Magenta)
    mask_sg1_and_sg2_v1 = is_sg1_v1 & is_sg2_v1 & ~is_ks_v1;   % Group 3: SG1 and SG2 (Green)
    mask_default_v1 = ~is_sg1_v1 & ~is_ks_v1;                  % Group 4: Default/Outside (Blue)
    
    % --- Classification for Vector 2 ---
    V2_Coeffs = vector2(:, 2);
    is_ks_v2 = ismember(V2_Coeffs, ks_coeff);
    is_sg1_v2 = ismember(V2_Coeffs, SG1_Coeffs);
    is_sg2_v2 = ismember(V2_Coeffs, SG2_Coeffs);
    
    % Mutually exclusive masks for V2
    mask_ks_v2 = is_ks_v2;
    mask_sg1_only_v2 = is_sg1_v2 & ~is_sg2_v2 & ~is_ks_v2;
    mask_sg1_and_sg2_v2 = is_sg1_v2 & is_sg2_v2 & ~is_ks_v2;
    mask_default_v2 = ~is_sg1_v2 & ~is_ks_v2;

    % Indices for plotting (X-axis, which is the sorted index)
    indices_v1 = 1:length(vector1);
    indices_v2 = 1:length(vector2);

    % --- 3. Plot Vector 1 Line and Markers (Left Y-Axis) ---
    yyaxis left;
    h_line1 = plot(vector1(:,1), '-', 'Color', COLOR_1, 'LineWidth', 1.5, 'DisplayName', plotID1);
    ylabel(sprintf('%s Value', plotID1));
    ax1 = gca;
    ax1.YColor = COLOR_1;
    
    % Plot V1 markers based on classification
    
    scatter(indices_v1(mask_ks_v1), vector1(mask_ks_v1, 1), 70, COLOR_KS, '^', 'filled', 'HandleVisibility', 'off');
    
    scatter(indices_v1(mask_sg1_only_v1), vector1(mask_sg1_only_v1, 1), 50, COLOR_SG1_ONLY, 'o', 'filled', 'HandleVisibility', 'off');
    
    scatter(indices_v1(mask_sg1_and_sg2_v1), vector1(mask_sg1_and_sg2_v1, 1), 50, COLOR_SG1_AND_SG2, 'o', 'filled', 'HandleVisibility', 'off');
    
    scatter(indices_v1(mask_default_v1), vector1(mask_default_v1, 1), 60, COLOR_1, 'x', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    
    yyaxis right;
    h_line2 = plot(vector2(:,1), '-', 'Color', COLOR_2, 'LineWidth', 1.5, 'DisplayName', plotID2);
    ylabel(sprintf('%s Value', plotID2));
    ax2 = gca;
    ax2.YColor = COLOR_2;
    
    % Plot V2 markers based on classification
    % Group 1: KS Coeff (Triangle - ^, Black)
    scatter(indices_v2(mask_ks_v2), vector2(mask_ks_v2, 1), 70, COLOR_KS, '^', 'filled', 'HandleVisibility', 'off');
    
    % Group 2: SG1 ONLY (Magenta Diamond - d)
    scatter(indices_v2(mask_sg1_only_v2), vector2(mask_sg1_only_v2, 1), 50, COLOR_SG1_ONLY, 'd', 'filled', 'HandleVisibility', 'off');
    
    % Group 3: SG1 AND SG2 (Green Diamond - d)
    scatter(indices_v2(mask_sg1_and_sg2_v2), vector2(mask_sg1_and_sg2_v2, 1), 50, COLOR_SG1_AND_SG2, 'd', 'filled', 'HandleVisibility', 'off');
    
    % Group 4: Default (Orange '+')
    scatter(indices_v2(mask_default_v2), vector2(mask_default_v2, 1), 60, COLOR_2, '+', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    
    yyaxis left; % Ensure X-axis lines are drawn on the correct axes
    
    % X-line for vector 1 cutoff
    h_xline1 = xline(cutoff1, '--', 'Color', COLOR_1, 'LineWidth', 1.5, 'Alpha', 0.8, 'DisplayName', [plotID1 ' Knee']);
    
    % X-line for vector 2 cutoff 
    h_xline2 = xline(cutoff2, ':', 'Color', COLOR_2, 'LineWidth', 1.5, 'Alpha', 0.8, 'DisplayName', [plotID2 ' Knee']);
    
    h_legend = [];
    
    % Line plots
    h_legend(end+1) = h_line1;
    h_legend(end+1) = h_line2;

    % Cutoff Lines
    h_legend(end+1) = h_xline1;
    h_legend(end+1) = h_xline2;
    
    % Classification Markers
    
    % KS Coeff (Triangle)
    h_legend(end+1) = scatter(NaN, NaN, 70, COLOR_KS, '^', 'filled', 'DisplayName', 'KS Coeff (Triangle)');

    % SG1 AND SG2 (Green Filled Marker)
    h_legend(end+1) = scatter(NaN, NaN, 50, COLOR_SG1_AND_SG2, 'o', 'filled', 'DisplayName', 'original GMM select');
    
    % SG1 ONLY (Magenta Filled Marker)
    h_legend(end+1) = scatter(NaN, NaN, 50, COLOR_SG1_ONLY, 'o', 'filled', 'DisplayName', 'removed by spike match');
    
    % Default/Outside (X marker)
    h_legend(end+1) = scatter(NaN, NaN, 60, COLOR_1, 'x', 'LineWidth', 1.5, 'DisplayName', 'Default/Outside');
    h_legend(end+1) = scatter(NaN, NaN, 60, COLOR_2, 'x', 'LineWidth', 1.5, 'DisplayName', 'Default/Outside');

    % Ensure handle array is valid and plot the legend
    h_legend = h_legend(isgraphics(h_legend));
    legend(h_legend, 'Location', 'bestoutside');
    
    hold off;
    
    if isnumeric(channelNum)
        channelNumStr = num2str(channelNum);
    else
        channelNumStr = channelNum;
    end
    
    if ~exist(folderName,'dir'), mkdir(folderName); end
    filename = fullfile(folderName, sprintf('ch%s_dual_classification_%s_vs_%s.pdf', channelNumStr, plotID1, plotID2));
    exportgraphics(fig, filename, 'Resolution', 300);
    close(fig);
end