function [threshQ2, threshQ4, select_gauss] = lineExclusion(medDist_vec, medDist_select, k_select, kDist_vec)
    % LINEEXCLUSION:
    % Draw Q2 and Q4 threshold lines limited to their quadrants,
    % match medDist and kDist by [coeff#, gauss#], and select
    % points above/right of each within their quadrants. Plots results.
    %
    % Output:
    %   threshQ2     - [slope, intercept] for Q2 line
    %   threshQ4     - [slope, intercept] for Q4 line
    %   select_gauss - [medDist, k_value, coeff#, gauss#] of selected points

    %% --- Compute limits
    maxKv = max(kDist_vec(:,1));
    minKv = min(kDist_vec(:,1));
    maxMedDist = max(medDist_vec(:,1));
    minMedDist = min(medDist_vec(:,1));
    
    % Intersection (shared vertex)
    x_inter = medDist_select(end,1);
    y_inter = k_select(end,1);

    %% --- Q2 line (intercept → minMedDist, maxKv)
    x1_Q2 = x_inter; y1_Q2 = y_inter;
    x2_Q2 = minMedDist; y2_Q2 = maxKv;
    slope_Q2 = (y2_Q2 - y1_Q2) / (x2_Q2 - x1_Q2);
    intercept_Q2 = y1_Q2 - slope_Q2 * x1_Q2;
    threshQ2 = [slope_Q2, intercept_Q2];

    %% --- Q4 line (intercept → maxMedDist, minKv)
    x1_Q4 = x_inter; y1_Q4 = y_inter;
    x2_Q4 = maxMedDist; y2_Q4 = minKv;
    slope_Q4 = (y2_Q4 - y1_Q4) / (x2_Q4 - x1_Q4);
    intercept_Q4 = y1_Q4 - slope_Q4 * x1_Q4;
    threshQ4 = [slope_Q4, intercept_Q4];

    %% --- Match medDist and kDist by [coeff#, gauss#]
    [commonPairs, idxMed, idxK] = intersect(medDist_vec(:,2:3), kDist_vec(:,2:3), 'rows');

    % Extract matched values
    x_vals = medDist_vec(idxMed, 1);
    y_vals = kDist_vec(idxK, 1);

    %% --- Compute line y-values
    y_line_Q2 = slope_Q2 * x_vals + intercept_Q2;
    y_line_Q4 = slope_Q4 * x_vals + intercept_Q4;

    %% --- Quadrant restriction
    in_Q2 = (x_vals <= x_inter) & (x_vals >= minMedDist);
    in_Q4 = (x_vals >= x_inter) & (x_vals <= maxMedDist);

    %% --- Selection logic (above/right of line within quadrant)
    select_Q2 = (y_vals > y_line_Q2) & in_Q2;
    select_Q4 = (y_vals > y_line_Q4) & in_Q4;

    select_mask = select_Q2 | select_Q4;

    %% --- Build full select_gauss table: [medDist, k_value, coeff#, gauss#]
    select_gauss = [x_vals(select_mask), y_vals(select_mask), commonPairs(select_mask,:)];

    % %% --- Plotting
    % figure; hold on; grid on; box on;
    % 
    % % Plot all matched points (gray)
    % scatter(x_vals, y_vals, 25, [0.6 0.6 0.6], 'filled');
    % 
    % % Highlight selected points (green)
    % scatter(x_vals(select_mask), y_vals(select_mask), 35, 'g', 'filled');
    % 
    % % Intersection
    % plot(x_inter, y_inter, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    % 
    % % Threshold lines (dashed black)
    % plot([x1_Q2, x2_Q2], [y1_Q2, y2_Q2], 'k--', 'LineWidth', 1.2);
    % plot([x1_Q4, x2_Q4], [y1_Q4, y2_Q4], 'k--', 'LineWidth', 1.2);
    % 
    % xlabel('Median Distance');
    % ylabel('k Distance');
    % title('Quadrant Threshold Lines and Selected Gaussians');
    % legend({'Matched points', 'Selected', 'Intersection', 'Threshold lines'}, 'Location', 'best');
    % hold off;
end
