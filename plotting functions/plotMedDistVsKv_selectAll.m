function plotMedDistVsKv_selectAll(medDist_vec, kDist_vec, medDist_select, kDist_select, folderName, channelNum)
% Plot Median Distance vs K-value with sequential boundary trimming
%
% This version:
% 1. Defines "Original" boundary using GLOBAL min/max of ALL points.
% 2. Calculates all distances from this GLOBAL boundary using a three-way split:
%    - Q2/Q3: Use distQ2
%    - Q4: Use distQ4
%    - Q1: Use min(distQ2, distQ4)
% 3. Sequentially trims the min/max points based on the distance cutoff to define new boundaries.
% 4. Selects points "above" the diagonal lines for each boundary.

if isnumeric(channelNum)
    channelNum = num2str(channelNum);
end
if ~exist(folderName,'dir')
    mkdir(folderName);
end

%% --- Match points by coeff# and gauss#
[commonPairs, idxMed, idxK] = intersect(medDist_vec(:,2:3), kDist_vec(:,2:3), 'rows');
medVals = medDist_vec(idxMed,1);
kVals   = kDist_vec(idxK,1);
coeffGauss = medDist_vec(idxMed,2:3);

%% --- Setup figure
figure('Visible','on'); hold on; grid on; box on;
scatter(medVals, kVals, 25, [0.7 0.7 0.7], 'filled', 'DisplayName', 'All Points'); % all points gray

% Add coeff# labels for all points
for i = 1:length(medVals)
    text(medVals(i), kVals(i)+0.7, num2str(coeffGauss(i,1)), ...
        'HorizontalAlignment','center', 'FontSize',7, 'FontWeight','bold', 'Color','k');
end

% Intersection
x_inter = medDist_select(end,1);
y_inter = kDist_select(end,1);
pt_inter = [x_inter, y_inter];
plot(x_inter, y_inter, 'ko', 'MarkerFaceColor','y', 'MarkerSize',8, 'DisplayName','Intersection');

%% --- Boundary setup
trimPerc = [10 5 2.5 1];
boundaryColors = {[0 0.5 0],[0 0 1],[0.5 0 0.5],[0 1 1]};
boundaryNames = {'pct10','pct5','pct2_5','pct1'};
boundaryDisplayNames = {'10%','5%','2.5%','1%'};

boundaryLines = struct();
allSelected = [];
h_legend = [];

%% --- 1. Define and Plot Original (0%) Boundary (Global Min/Max)
pt_Q2_orig_vertex = [min(medVals), max(kVals)];
pt_Q4_orig_vertex = [max(medVals), min(kVals)];

% Plot the original boundary lines
h_orig_q2 = plot([pt_inter(1), pt_Q2_orig_vertex(1)], [pt_inter(2), pt_Q2_orig_vertex(2)], '-', 'Color', [0 0 0], 'LineWidth',1.5, 'DisplayName', 'Original Boundary');
plot([pt_inter(1), pt_Q4_orig_vertex(1)], [pt_inter(2), pt_Q4_orig_vertex(2)], '-', 'Color', [0 0 0], 'LineWidth',1.5, 'HandleVisibility','off');
h_legend(end+1) = h_orig_q2;

% Store in struct
boundaryLines.orig = struct('Q2',[pt_inter; pt_Q2_orig_vertex],'Q4',[pt_inter; pt_Q4_orig_vertex]);

%% --- 2. Calculate Distances to ORIGINAL Boundary (Q1 min-distance rule applied here)
% Get line equations for original boundary: y = m*x + c
% Q2 Line
m1 = (pt_Q2_orig_vertex(2) - pt_inter(2)) / (pt_Q2_orig_vertex(1) - pt_inter(1));
c1 = pt_inter(2) - m1 * pt_inter(1);
if isinf(m1) % Handle vertical line case
    m1 = realmax; c1 = 0; 
end

% Q4 Line
m2 = (pt_Q4_orig_vertex(2) - pt_inter(2)) / (pt_Q4_orig_vertex(1) - pt_inter(1));
c2 = pt_inter(2) - m2 * pt_inter(1);
if isinf(m2) % Handle vertical line case
    m2 = realmax; c2 = 0;
end

% Calculate perpendicular distance of ALL points to the ORIGINAL boundary
distQ2 = abs(m1 * medVals - kVals + c1) / sqrt(m1^2 + 1);
distQ4 = abs(m2 * medVals - kVals + c2) / sqrt(m2^2 + 1);

% Apply the NEW distance logic:
dist_to_orig_boundary = zeros(size(medVals));

% Left of intersection (Q2/Q3): Use distQ2
idx_left = medVals <= x_inter;
dist_to_orig_boundary(idx_left) = distQ2(idx_left);

% Right of intersection (Q1/Q4):
idx_right = medVals > x_inter;

% Q1 Points: medVals > x_inter AND kVals > y_inter
idx_Q1 = idx_right & (kVals > y_inter);
dist_to_orig_boundary(idx_Q1) = min(distQ2(idx_Q1), distQ4(idx_Q1));

% Q4 Points: medVals > x_inter AND kVals <= y_inter
idx_Q4 = idx_right & (kVals <= y_inter);
dist_to_orig_boundary(idx_Q4) = distQ4(idx_Q4); % Use Q4 distance

%% --- 3. Calculate and Plot Trimmed Boundaries
for t = 1:length(trimPerc)
    pct = trimPerc(t);
    
    % 1. Determine "kept" points for boundary creation
    cutoff = prctile(dist_to_orig_boundary, 100 - pct);
    keepMask = (dist_to_orig_boundary <= cutoff);
    
    med_keep = medVals(keepMask);
    k_keep   = kVals(keepMask);
    
    % 2. Recompute boundary vertices based on the GLOBAL min/max of 
    % the *entire* set of kept points.
    if ~isempty(med_keep)
        Q2_vertex_new = [min(med_keep), max(k_keep)];
        Q4_vertex_new = [max(med_keep), min(k_keep)];
    else
        Q2_vertex_new = pt_inter; % Fallback
        Q4_vertex_new = pt_inter; % Fallback
    end
    
    Q2_pt = [pt_inter; Q2_vertex_new];
    Q4_pt = [pt_inter; Q4_vertex_new];

    % 3. Determine selected points based on NEW boundary LINES (same as before)
    m_Q2 = (Q2_pt(2,2) - Q2_pt(1,2)) / (Q2_pt(2,1) - Q2_pt(1,1));
    c_Q2 = Q2_pt(1,2) - m_Q2 * Q2_pt(1,1);
    if abs(Q2_pt(2,1) - Q2_pt(1,1)) < 1e-9 
        m_Q2 = Inf; 
        k_thresh_Q2 = Inf(size(medVals));
    else
        k_thresh_Q2 = m_Q2 * medVals + c_Q2;
    end

    m_Q4 = (Q4_pt(2,2) - Q4_pt(1,2)) / (Q4_pt(2,1) - Q4_pt(1,1));
    c_Q4 = Q4_pt(1,2) - m_Q4 * Q4_pt(1,1);
    if abs(Q4_pt(2,1) - Q4_pt(1,1)) < 1e-9
        m_Q4 = Inf;
        k_thresh_Q4 = Inf(size(medVals));
    else
        k_thresh_Q4 = m_Q4 * medVals + c_Q4;
    end

    % Select from ALL original points
    selMask = ((medVals <= x_inter) & (kVals >= k_thresh_Q2)) | ...
              ((medVals > x_inter)  & (kVals >= k_thresh_Q4));
          
    selected = [medVals(selMask), kVals(selMask), coeffGauss(selMask,:)];
    allSelected = [allSelected; selected(:,3:4)];
    
    boundaryLines.(boundaryNames{t}) = struct('Q2',Q2_pt,'Q4',Q4_pt,'selected',selected);
    
    % --- Plot boundary lines
    h = plot([Q2_pt(1,1), Q2_pt(2,1)], [pt_inter(2), Q2_pt(2,2)], '-', 'Color', boundaryColors{t}, 'LineWidth',1.5);
    plot([pt_inter(1), Q4_pt(2,1)], [pt_inter(2), Q4_pt(2,2)], '-', 'Color', boundaryColors{t}, 'LineWidth',1.5, 'HandleVisibility','off');
    
    % Add to legend
    h.DisplayName = sprintf('Boundary %s', boundaryDisplayNames{t});
    h_legend(end+1) = h;
    
    % --- Plot selected points
    scatter(selected(:,1), selected(:,2), 50, boundaryColors{t}, 'filled', 'HandleVisibility','off');
end

% --- 4. Select and plot points for the ORIGINAL boundary ---
m_Q2_orig = m1; c_Q2_orig = c1;
m_Q4_orig = m2; c_Q4_orig = c2;

k_thresh_Q2_orig = m_Q2_orig * medVals + c_Q2_orig;
k_thresh_Q4_orig = m_Q4_orig * medVals + c_Q4_orig;

selMask_orig = ((medVals <= x_inter) & (kVals >= k_thresh_Q2_orig)) | ...
               ((medVals > x_inter)  & (kVals >= k_thresh_Q4_orig));

selected_orig = [medVals(selMask_orig), kVals(selMask_orig), coeffGauss(selMask_orig,:)];
allSelected = [allSelected; selected_orig(:,3:4)];
boundaryLines.orig.selected = selected_orig;

% Plot selected points for original boundary (black)
scatter(selected_orig(:,1), selected_orig(:,2), 50, [0 0 0], 'filled', 'HandleVisibility','off');

%% --- Plot remaining points
% Make sure we don't re-plot selected points
allSelected = unique(allSelected, 'rows');
plotted_mask = ismember(coeffGauss, allSelected, 'rows');

for k = 1:length(medVals)
    if plotted_mask(k)
        continue % already plotted as selected
    end
    
    coeff_num = coeffGauss(k,1);
    
    % same coeff, different gauss
    if ismember(coeff_num, allSelected(:,1))
        h_same = scatter(medVals(k), kVals(k), 35, [0.5 0.5 1], 'filled');
    else
        h_unrelated = scatter(medVals(k), kVals(k), 35, [1 0 0], 'x', 'LineWidth',1.5);
    end
end

%% --- Legend
if exist('h_same', 'var')
    h_same.DisplayName = 'Same coeff';
    h_legend(end+1) = h_same(1);
end
if exist('h_unrelated', 'var')
    h_unrelated.DisplayName = 'Unrelated';
    h_legend(end+1) = h_unrelated(1);
end

legend(h_legend,'Location','bestoutside');
xlabel('Median Distance');
ylabel('K-value');
title('Median Distance vs K-value with Sequential Boundaries');

%% --- Save figure
filename = fullfile(folderName, sprintf('ch%s_medD_v_kv_boundaries.png', channelNum));
exportgraphics(gcf, filename, 'Resolution',300);
hold off;
end