function pltDualKneePlot(len1, vector1, len2, vector2, plotID1, plotID2, channelNum, folderName)
% PLTDUALKNEEPLOT plots two vectors on a dual Y-axis figure, applying a
% "knee" selection cutoff to both.
%
% Inputs:
%   len1:       Length of the 'selected' portion of vector1 (from the end).
%   vector1:    Full sorted vector 1: [value, coeff#, gauss#]
%   len2:       Length of the 'selected' portion of vector2 (from the end).
%   vector2:    Full sorted vector 2: [value, coeff#, gauss#]
%   plotID1:    String identifier for vector 1 (e.g., 'kv').
%   plotID2:    String identifier for vector 2 (e.g., 'medDist').
%   channelNum: Channel number for filename.
%   folderName: Folder to save the output file.

    % Define colors
    COLOR_1 = [0.0 0.447 0.741]; % MATLAB Blue for Left Axis
    COLOR_2 = [0.85 0.325 0.098]; % MATLAB Red/Orange for Right Axis

    % --- 1. Setup Figure and Axes ---
    fig = figure('Visible', 'on');
    hold on;
    grid on;
    
    title(sprintf('%s vs %s per Coefficient (Channel %s)', plotID1, plotID2, num2str(channelNum)));
    xlabel('Coefficient Index (Sorted)');
    
    % Get total length for x-axis
    max_len = max(length(vector1), length(vector2));
    
    % Calculate cutoff points
    cutoff1 = length(vector1) - len1;
    cutoff2 = length(vector2) - len2;

    % --- 2. Plot Vector 1 (Left Y-Axis) ---
    yyaxis left;
    h_line1 = plot(vector1(:,1), '-', 'Color', COLOR_1, 'LineWidth', 1.5, 'DisplayName', plotID1);
    ylabel(sprintf('%s Value', plotID1));
    ax1 = gca;
    ax1.YColor = COLOR_1;
    
    % Plot selected/not-selected markers for vector 1
    h_select1 = [];
    h_notselect1 = [];
    
    % Use logical indexing for cleaner marker plotting
    indices = 1:length(vector1);
    
    % Selected points (right side of cutoff)
    selected_mask1 = indices > cutoff1;
    h_select1 = scatter(indices(selected_mask1), vector1(selected_mask1, 1), 50, COLOR_1, 'o', 'filled', 'HandleVisibility', 'off');
    
    % Not selected points (left side of cutoff)
    not_selected_mask1 = indices <= cutoff1;
    h_notselect1 = scatter(indices(not_selected_mask1), vector1(not_selected_mask1, 1), 60, COLOR_1, 'x', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % --- 3. Plot Vector 2 (Right Y-Axis) ---
    yyaxis right;
    h_line2 = plot(vector2(:,1), '-', 'Color', COLOR_2, 'LineWidth', 1.5, 'DisplayName', plotID2);
    ylabel(sprintf('%s Value', plotID2));
    ax2 = gca;
    ax2.YColor = COLOR_2;
    
    % Plot selected/not-selected markers for vector 2
    h_select2 = [];
    h_notselect2 = [];

    % Use logical indexing for cleaner marker plotting
    indices = 1:length(vector2);
    
    % Selected points (right side of cutoff)
    selected_mask2 = indices > cutoff2;
    h_select2 = scatter(indices(selected_mask2), vector2(selected_mask2, 1), 50, COLOR_2, 'd', 'filled', 'HandleVisibility', 'off'); % Use diamond 'd' for distinction
    
    % Not selected points (left side of cutoff)
    not_selected_mask2 = indices <= cutoff2;
    h_notselect2 = scatter(indices(not_selected_mask2), vector2(not_selected_mask2, 1), 60, COLOR_2, '+', 'LineWidth', 1.5, 'HandleVisibility', 'off'); % Use plus '+' for distinction

    % --- 4. Add Reference Lines (X-axis must be consistent) ---
    % Switch back to left axis for shared elements like xline
    yyaxis left;
    
    % X-line for vector 1 cutoff
    h_xline1 = xline(cutoff1, '--', 'Color', COLOR_1, 'LineWidth', 1.5, 'Alpha', 0.8, 'DisplayName', [plotID1 ' Cutoff']);
    
    % X-line for vector 2 cutoff (Use a different style or color for visual separation)
    h_xline2 = xline(cutoff2, ':', 'Color', COLOR_2, 'LineWidth', 1.5, 'Alpha', 0.8, 'DisplayName', [plotID2 ' Cutoff']);

    % --- 5. Build Legend ---
    h_legend = [];
    
    % Group 1 (Line 1, Selected 1, Not Selected 1, Cutoff 1)
    h_legend(end+1) = h_line1;
    h_legend(end+1) = scatter(NaN, NaN, 50, COLOR_1, 'o', 'filled', 'DisplayName', [plotID1 ' Selected']);
    h_legend(end+1) = scatter(NaN, NaN, 60, COLOR_1, 'x', 'LineWidth', 1.5, 'DisplayName', [plotID1 ' Not Selected']);
    h_legend(end+1) = h_xline1;
    
    % Group 2 (Line 2, Selected 2, Not Selected 2, Cutoff 2)
    h_legend(end+1) = h_line2;
    h_legend(end+1) = scatter(NaN, NaN, 50, COLOR_2, 'd', 'filled', 'DisplayName', [plotID2 ' Selected']);
    h_legend(end+1) = scatter(NaN, NaN, 60, COLOR_2, '+', 'LineWidth', 1.5, 'DisplayName', [plotID2 ' Not Selected']);
    h_legend(end+1) = h_xline2;
    
    % Ensure handle array is valid and plot the legend
    h_legend = h_legend(isgraphics(h_legend));
    legend(h_legend, 'Location', 'bestoutside');
    
    hold off;

    % --- 6. Save Figure ---
    if isnumeric(channelNum)
        channelNumStr = num2str(channelNum);
    else
        channelNumStr = channelNum;
    end
    
    if ~exist(folderName,'dir'), mkdir(folderName); end

    filename = fullfile(folderName, sprintf('ch%s_dual_knee_%s_vs_%s.png', channelNumStr, plotID1, plotID2));
    exportgraphics(fig, filename, 'Resolution', 300);
    close(fig);
end