function [select_gauss_orig, select_gauss_1pct, select_gauss_2_5pct, select_gauss_5pct, select_gauss_10pct, select_all] = lineExclusion(medDist_vec, medDist_select, k_select, kDist_vec, folderName, channelNum)
% LINEEXCLUSION (Final: Updated with Quadrant-Specific Trimming Logic and Two Optional Plots)
% Builds Q2/Q4 threshold lines, iteratively trims based on quadrant-specific distance,
% and returns selected Gaussians while optionally saving two figures.

% --- Default/Optional Input Handling ---
if nargin < 5, folderName = ''; end 
if nargin < 6, channelNum = ''; end 

plot_and_save = ~isempty(folderName) && ~isempty(channelNum); 

%% --- INITIALIZATION & DATA MATCHING ---
maxKv = max(kDist_vec(:,1));
minKv = min(kDist_vec(:,1));
maxMedDist = max(medDist_vec(:,1));
minMedDist = min(medDist_vec(:,1));
x_inter = medDist_select(end,1);
y_inter = k_select(end,1);
pt_inter = [x_inter, y_inter]; 

[commonPairs, idxMed, idxK] = intersect(medDist_vec(:,2:3), kDist_vec(:,2:3), 'rows');
x_vals = medDist_vec(idxMed, 1);
y_vals = kDist_vec(idxK, 1);
gauss_ids = commonPairs; 
coeffGauss = medDist_vec(idxMed,2:3); 

%% --- ORIGINAL BOUNDARY & DISTANCE CALCULATION (FOR TRIMMING) ---
% 1. Original Boundaries (0% trim)
slope_Q2 = (maxKv - y_inter) / (minMedDist - x_inter);
slope_Q4 = (minKv - y_inter) / (maxMedDist - x_inter);
int_Q2 = y_inter - slope_Q2 * x_inter;
int_Q4 = y_inter - slope_Q4 * x_inter;

% 2. Get line parameters for distance calculation (Handling vertical lines)
m1 = slope_Q2; c1 = int_Q2;
m2 = slope_Q4; c2 = int_Q4;
if abs(minMedDist - x_inter) < 1e-9, m1 = realmax; c1 = 0; end
if abs(maxMedDist - x_inter) < 1e-9, m2 = realmax; c2 = 0; end

% 3. Calculate perpendicular distance of ALL points to the ORIGINAL boundary
distQ2 = abs(m1 * x_vals - y_vals + c1) / sqrt(m1^2 + 1);
distQ4 = abs(m2 * x_vals - y_vals + c2) / sqrt(m2^2 + 1);

% 4. APPLY QUADRANT-SPECIFIC DISTANCE LOGIC FOR TRIMMING
dist_to_boundary = zeros(size(x_vals));
idx_right = x_vals > x_inter;

% Q1 Points: min(distQ2, distQ4)
idx_Q1 = idx_right & (y_vals > y_inter);
dist_to_boundary(idx_Q1) = min(distQ2(idx_Q1), distQ4(idx_Q1));

% Q4 Points: distQ4
idx_Q4 = idx_right & (y_vals <= y_inter);
dist_to_boundary(idx_Q4) = distQ4(idx_Q4); 

% Q2/Q3 Points: distQ2
idx_left = x_vals <= x_inter;
dist_to_boundary(idx_left) = distQ2(idx_left);

%% --- ORIGINAL SELECTION (0% trim) ---
y_line_Q2 = slope_Q2*x_vals + int_Q2;
y_line_Q4 = slope_Q4*x_vals + int_Q4;
in_Q2 = (x_vals <= x_inter) & (x_vals >= minMedDist);
in_Q4 = (x_vals >= x_inter) & (x_vals <= maxMedDist);

% Use inclusive inequality (>=) for selection
select_mask_orig = ((y_vals >= y_line_Q2) & in_Q2) | ((y_vals >= y_line_Q4) & in_Q4);
select_gauss_orig = [x_vals(select_mask_orig), y_vals(select_mask_orig), gauss_ids(select_mask_orig,:)];

select_all.original = struct('threshQ2',[slope_Q2,int_Q2],...
                             'threshQ4',[slope_Q4,int_Q4],...
                             'select_gauss',select_gauss_orig,...
                             'removed_points',0,...
                             'keep_mask', true(size(x_vals)));

%% --- SEQUENTIAL TRIMMING AND RE-SELECTION ---
trim_list = [1, 2.5, 5, 10];
colors = {[0 0 1],[0.5 0 0.5],[0 1 1],[0 0.5 0]};
field_names = {'trim1_0pct','trim2_5pct','trim5_0pct','trim10_0pct'};
select_sets = cell(1,numel(trim_list));
allSelectedGauss = [select_gauss_orig(:,3:4)]; 

for ii = 1:numel(trim_list)
    trimPct = trim_list(ii);
    
    % 1. Determine points to keep for boundary re-fitting
    cutoff = prctile(dist_to_boundary, 100 - trimPct);
    keep_mask = dist_to_boundary <= cutoff;
    removed_mask = ~keep_mask;
    x_keep = x_vals(keep_mask);
    y_keep = y_vals(keep_mask);
    
    % 2. Refit thresholds 
    slope_Q2r = (max(y_keep) - y_inter) / (min(x_keep) - x_inter);
    slope_Q4r = (min(y_keep) - y_inter) / (max(x_keep) - x_inter);
    int_Q2r = y_inter - slope_Q2r*x_inter;
    int_Q4r = y_inter - slope_Q4r*x_inter;
    
    % 3. Re-select (using ALL points against the new boundaries, inclusive >=)
    y_line_Q2r = slope_Q2r*x_vals + int_Q2r;
    y_line_Q4r = slope_Q4r*x_vals + int_Q4r;
    in_Q2r = (x_vals <= x_inter) & (x_vals >= minMedDist);
    in_Q4r = (x_vals >= x_inter) & (x_vals <= maxMedDist);
    select_mask_r = ((y_vals >= y_line_Q2r) & in_Q2r) | ((y_vals >= y_line_Q4r) & in_Q4r);
    select_gauss_refined = [x_vals(select_mask_r), y_vals(select_mask_r), gauss_ids(select_mask_r,:)];
    
    select_all.(field_names{ii}) = struct('threshQ2',[slope_Q2r,int_Q2r],...
                                          'threshQ4',[slope_Q4r,int_Q4r],...
                                          'select_gauss',select_gauss_refined,...
                                          'removed_points',sum(removed_mask),...
                                          'keep_mask', keep_mask);
    select_sets{ii} = select_gauss_refined;
    allSelectedGauss = [allSelectedGauss; select_gauss_refined(:,3:4)];
end

% Assign individual outputs
select_gauss_1pct   = select_sets{1};
select_gauss_2_5pct = select_sets{2};
select_gauss_5pct   = select_sets{3};
select_gauss_10pct  = select_sets{4};

%% --- OPTIONAL PLOTTING AND SAVING TWO FIGURES ---
if plot_and_save
    
    if isnumeric(channelNum), channelNumStr = num2str(channelNum); else channelNumStr = channelNum; end
    if ~exist(folderName,'dir'), mkdir(folderName); end
    
    x_left = linspace(minMedDist, x_inter, 300);
    x_right = linspace(x_inter, maxMedDist, 300);
    
    % ---------------------------------------------------------------
    % FIGURE 1: Boundary & Selection Plot (Detailed Labeling)
    % ---------------------------------------------------------------
    figure('Visible','on'); hold on; grid on; box on;
    h_legend = [];
    
    % Plot ALL points (gray background)
    h_all = scatter(x_vals, y_vals, 25, [0.7 0.7 0.7], 'filled', 'DisplayName', 'All Points'); 
    h_legend(end+1) = h_all;
    
    % Add coeff# labels for all points
    for i = 1:length(x_vals)
        text(x_vals(i), y_vals(i)+0.7, num2str(coeffGauss(i,1)), ...
            'HorizontalAlignment','center', 'FontSize',7, 'FontWeight','bold', 'Color','k');
    end
    
    % Plot Intersection point (Yellow circle)
    h_inter = plot(x_inter, y_inter, 'ko', 'MarkerFaceColor','y', 'MarkerSize',8, 'DisplayName','Intersection');
    h_legend(end+1) = h_inter;

    % Plot Trimmed Boundaries (1%, 2.5%, 5%, 10%)
    boundaryDisplayNames = {'1%','2.5%','5%','10%'};
    boundaryColors = {[0 0 1],[0.5 0 0.5],[0 1 1],[0 0.5 0]};

    for ii = 1:numel(trim_list)
        trim_struct = select_all.(field_names{ii});
        sQ2 = trim_struct.threshQ2(1); iQ2 = trim_struct.threshQ2(2);
        sQ4 = trim_struct.threshQ4(1); iQ4 = trim_struct.threshQ4(2);
        color = colors{ii};
    
        % Use the points that were actually used to calculate this boundary
        x_keep = x_vals(trim_struct.keep_mask);
        y_keep = y_vals(trim_struct.keep_mask);
    
        % Compute intersection for this trimmed boundary
        x_inter_trim = x_keep(end);  % or wherever your intersection logic sets it
        y_inter_trim = y_keep(end);
    
        % Q2 boundary segment (left of intercept)
        xQ2_pts = x_keep(x_keep <= x_inter_trim);  % points used to fit Q2
        if ~isempty(xQ2_pts)
            x_start = min(xQ2_pts);      % leftmost boundary point
            x_end   = x_inter_trim;      % intersection
            y_start = slope_Q2r*x_start + int_Q2r;
            y_end   = slope_Q2r*x_end   + int_Q2r;
            plot([x_start, x_end], [y_start, y_end], '-', 'Color', color, 'LineWidth',1.5);
        end
        
        % Q4 boundary segment (right of intercept)
        xQ4_pts = x_keep(x_keep >= x_inter_trim);  % points used to fit Q4
        if ~isempty(xQ4_pts)
            x_start = x_inter_trim;       % intersection
            x_end   = max(xQ4_pts);       % rightmost boundary point
            y_start = slope_Q4r*x_start + int_Q4r;
            y_end   = slope_Q4r*x_end   + int_Q4r;
            plot([x_start, x_end], [y_start, y_end], '-', 'Color', color, 'LineWidth',1.5);
        end

    
        scatter(trim_struct.select_gauss(:,1), trim_struct.select_gauss(:,2), 50, color, 'filled', 'HandleVisibility','off');
    end



    % Plot Original Boundary (0%) and selected points (Black)
    sQ2_orig=select_all.original.threshQ2(1); iQ2_orig=select_all.original.threshQ2(2);
    sQ4_orig=select_all.original.threshQ4(1); iQ4_orig=select_all.original.threshQ4(2);
    
    xQ2_orig = linspace(min(x_vals), x_inter, 300);
    xQ4_orig = linspace(x_inter, max(x_vals), 300);
    
    h_orig_q2 = plot(xQ2_orig, sQ2_orig*xQ2_orig+iQ2_orig, '-', 'Color', [0 0 0], 'LineWidth',1.5, 'DisplayName', 'Original Boundary');
    plot(xQ4_orig, sQ4_orig*xQ4_orig+iQ4_orig, '-', 'Color', [0 0 0], 'LineWidth',1.5, 'HandleVisibility','off');

    h_legend(end+1) = h_orig_q2;
    selected_orig_pts = select_all.original.select_gauss;
    scatter(selected_orig_pts(:,1), selected_orig_pts(:,2), 50, [0 0 0], 'filled', 'HandleVisibility','off');
    
    
    % Plot remaining points (Same Coeff/Unrelated)
    allSelectedGauss = unique(allSelectedGauss, 'rows');
    plotted_mask = ismember(coeffGauss, allSelectedGauss, 'rows');
    h_same = []; h_unrelated = [];
    
    for k = 1:length(x_vals)
        if plotted_mask(k), continue; end
        coeff_num = coeffGauss(k,1);
        if ismember(coeff_num, allSelectedGauss(:,1))
            h_same = scatter(x_vals(k), y_vals(k), 35, [0.5 0.5 1], 'filled');
        else
            h_unrelated = scatter(x_vals(k), y_vals(k), 35, [1 0 0], 'x', 'LineWidth',1.5);
        end
    end
    
    % Final Legend Cleanup
    if exist('h_same', 'var') && ~isempty(h_same), h_same.DisplayName = 'Same coeff'; h_legend(end+1) = h_same(1); end
    if exist('h_unrelated', 'var') && ~isempty(h_unrelated), h_unrelated.DisplayName = 'Unrelated'; h_legend(end+1) = h_unrelated(1); end
    
    legend(h_legend,'Location','bestoutside');
    % --- Add text under legend showing summary counts ---
    all_selected = select_all.original.select_gauss;
    num_gauss = size(all_selected, 1);
    num_unique_coeff = numel(unique(all_selected(:,3)));
    
    % Create annotation textbox below `
    txt_summary = sprintf('Selected Gaussians: %d\nUnique Coefficients: %d', num_gauss, num_unique_coeff);
    annotation('textbox',[0.83, 0.15, 0.15, 0.07], 'String',txt_summary, ...
    'FitBoxToText','on', 'EdgeColor','none', 'FontSize',9, 'VerticalAlignment','top');

    xlabel('Median Distance');
    ylabel('K-value');
    title('Median Distance vs K-value with Sequential Boundaries');
    
    % Save figure 1
    filename1 = fullfile(folderName, sprintf('ch%s_medD_v_kv_boundaries.png', channelNumStr));
    exportgraphics(gcf, filename1, 'Resolution',300);
    close(gcf); % Close the current figure

    % ---------------------------------------------------------------
    % FIGURE 2: Removal Plot (Shows points removed for trimming)
    % ---------------------------------------------------------------
    figure('Visible','on'); hold on; grid on; box on;

    % Use the 10% mask to define removed points
    removed_10pct_mask = ~select_all.trim10_0pct.keep_mask;
    
    % Plot all points that were kept (gray)
    h_kept = scatter(x_vals(~removed_10pct_mask), y_vals(~removed_10pct_mask), 18, [0.7 0.7 0.7], 'filled', 'DisplayName','Kept Points (10% Trim)');
    
    % Plot removed points (red 'x')
    h_removed = scatter(x_vals(removed_10pct_mask), y_vals(removed_10pct_mask), 30, [1 0 0], 'x', 'LineWidth', 1.5, 'DisplayName','Removed Points (10% Trim)');

    % Plot Boundaries (Original and Final 10% trim)
    plot(x_left, sQ2_orig*x_left+iQ2_orig, 'k--', 'LineWidth',1.4, 'DisplayName','Original Boundary');
    plot(x_right, sQ4_orig*x_right+iQ4_orig, 'k--', 'LineWidth',1.4, 'HandleVisibility','off');

    sQ2_10=select_all.trim10_0pct.threshQ2(1); iQ2_10=select_all.trim10_0pct.threshQ2(2);
    sQ4_10=select_all.trim10_0pct.threshQ4(1); iQ4_10=select_all.trim10_0pct.threshQ4(2);
    plot(x_left, sQ2_10*x_left+iQ2_10, 'g-', 'LineWidth',1.5, 'DisplayName','10% Trimmed Boundary');
    plot(x_right, sQ4_10*x_right+iQ4_10, 'g-', 'LineWidth',1.5, 'HandleVisibility','off');

    plot(x_inter, y_inter, 'ko', 'MarkerFaceColor','y', 'MarkerSize',8, 'DisplayName','Intersection','HandleVisibility','off');

    xlabel('Median Distance');
    ylabel('K-value');
    title(sprintf('Points Removed by 10%% Distance Cutoff (Channel %s)', channelNumStr));
    legend([h_kept, h_removed],'Location','bestoutside');

    % Save figure 2
    filename2 = fullfile(folderName, sprintf('ch%s_lineExclusion_removed_points.png', channelNumStr));
    exportgraphics(gcf, filename2, 'Resolution',300);
    close(gcf);
end

%% --- Print summary
    fprintf('\n=== LineExclusion Summary ===\n');
    fprintf('Original: %d selected\n', size(select_gauss_orig,1));
    fprintf('Trim 1%%: %d selected, %d removed\n', size(select_gauss_1pct,1), select_all.trim1_0pct.removed_points);
    fprintf('Trim 2.5%%: %d selected, %d removed\n', size(select_gauss_2_5pct,1), select_all.trim2_5pct.removed_points);
    fprintf('Trim 5%%: %d selected, %d removed\n', size(select_gauss_5pct,1), select_all.trim5_0pct.removed_points);
    fprintf('Trim 10%%: %d selected, %d removed\n', size(select_gauss_10pct,1), select_all.trim10_0pct.removed_points);
    fprintf('=============================\n\n');
end