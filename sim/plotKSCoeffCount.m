function plotKSCoeffCount(ks_coeff_all, select_all, unitNumbers, plotTitle, graphNum)
% PLOTKSCOEFFCOUNT plots number of coefficients per unit with shading, jitter,
% trend lines, R^2, and Spearman correlation.
%
% Inputs:
%   ks_coeff_all : cell array of ks coefficients per unit
%   select_all   : cell array of select_all coefficients per unit (GMM)
%   unitNumbers  : vector of unit numbers (X-axis)
%   plotTitle    : string, title for the plot
%   graphNum     : integer (1-5) to define suffix for saving

numUnits = numel(unitNumbers);
figure('Position', [100, 100, 1000, 750]); % Make figure larger
hold on;

% Ensure unitNumbers is a column vector for calculations
unitNumbers = unitNumbers(:); 
ks_counts  = cellfun(@numel, ks_coeff_all);
ks_counts = ks_counts(:); % Ensure column vector
sel_counts = cellfun(@(x) numel(unique(x)), select_all);  % unique counts (GMM)
sel_counts = sel_counts(:); % Ensure column vector

% Determine max count for axis limits
yMax = ceil(max([ks_counts; sel_counts]) * 1.2); % take max, scale 20%
if yMax == 0 % Handle case where all counts are 0
    yMax = 10;
end

for i = 1:numUnits
    unit = unitNumbers(i);
    if mod(unit,2) == 0
        xPatch = [unit-0.5, unit+0.5, unit+0.5, unit-0.5];
        yPatch = [0, 0, yMax, yMax]; % Use dynamic yMax
        patch(xPatch, yPatch, [0.95 0.95 0.95], 'EdgeColor','none'); % light gray
    end
end

hKS = [];
hSel = [];
jitterAmp = 0.5;  % +/- 0.1
ks_color = [0.2 0.6 0.8];  % Blue
sel_color = [0.8 0.4 0.2]; % Red-orange

markerSize_new = 6; % Changed MarkerSize here

for i = 1:numUnits
    x_jitter = unitNumbers(i) + (rand()*jitterAmp - jitterAmp/2);
    
    % KS coefficients
    ks_count = ks_counts(i);
    if isempty(hKS)
        hKS = plot(x_jitter, ks_count, 'o', 'MarkerFaceColor', ks_color, 'MarkerEdgeColor','k', 'MarkerSize', markerSize_new);
    else
        plot(x_jitter, ks_count, 'o', 'MarkerFaceColor', ks_color, 'MarkerEdgeColor','k', 'MarkerSize', markerSize_new);
    end
    
    % select_all coefficients (GMM)
    sel_count = sel_counts(i);
    if isempty(hSel)
        hSel = plot(x_jitter, sel_count, 's', 'MarkerFaceColor', sel_color, 'MarkerEdgeColor','k', 'MarkerSize', markerSize_new);
    else
        plot(x_jitter, sel_count, 's', 'MarkerFaceColor', sel_color, 'MarkerEdgeColor','k', 'MarkerSize', markerSize_new);
    end
end


% KS Data (Blue)
p_ks = polyfit(unitNumbers, ks_counts, 1);       % Linear fit
y_fit_ks = polyval(p_ks, unitNumbers);           % Y-values for trend line
[r_ks, ~] = corr(unitNumbers, ks_counts, 'type', 'Pearson');
R2_ks = r_ks^2;                                  % R-squared
[rho_ks, pval_ks] = corr(unitNumbers, ks_counts, 'type', 'Spearman');
% Plot KS trend line
hKS_fit = plot(unitNumbers, y_fit_ks, '-', 'Color', ks_color, 'LineWidth', 4);

% GMM Data (Red-orange)
p_sel = polyfit(unitNumbers, sel_counts, 1);     % Linear fit
y_fit_sel = polyval(p_sel, unitNumbers);         % Y-values for trend line
[r_sel, ~] = corr(unitNumbers, sel_counts, 'type', 'Pearson');
R2_sel = r_sel^2;                                % R-squared
[rho_sel, pval_sel] = corr(unitNumbers, sel_counts, 'type', 'Spearman');
% Plot GMM trend line
hSel_fit = plot(unitNumbers, y_fit_sel, '--', 'Color', sel_color, 'LineWidth', 4);


ax = gca; % Get current axes
ax.FontSize = 18; % Set tick/label font size

% Labels
xlabel('Unit Number', 'FontSize', 15, 'FontWeight', 'bold');
ylabel('Number of Coefficients', 'FontSize', 15, 'FontWeight', 'bold');
title(plotTitle, 'FontSize', 18, 'FontWeight', 'bold');

% Limits and Grid
xlim([min(unitNumbers)-0.5, max(unitNumbers)+0.5]);
ylim([0 yMax]);
grid on;
box on; % Add a box around the plot
ax.Layer = 'top'; % Ensure grid/ticks are on top of the patch

% Create legend strings with stats
% Using tex interpreter for rho (ρ) and italic p
ks_fit_str = sprintf('KS Trend (R^2=%.2f, {\\it\\rho}=%.2f, {\\itp}=%.3g)', R2_ks, rho_ks, pval_ks);
gmm_fit_str = sprintf('GMM Trend (R^2=%.2f, {\\it\\rho}=%.2f, {\\itp}=%.3g)', R2_sel, rho_sel, pval_sel);

% Create legend
hL = legend([hKS, hSel, hKS_fit, hSel_fit], ...
    {'KS Coefficients', 'GMM Coefficients', ks_fit_str, gmm_fit_str}, ...
    'Location', 'northwest', 'FontSize', 13);
hL.Interpreter = 'tex'; % Enable tex for formatting

folderSave = fullfile('results','counts');
if ~exist(folderSave,'dir'), mkdir(folderSave); end
switch graphNum
    case 1, suffix = '_0pct';
    case 2, suffix = '_1pct';
    case 3, suffix = '_2p5pct';
    case 4, suffix = '_5pct';
    case 5, suffix = '_10pct';
    otherwise, suffix = 'spikeMatch';
end

saveas(gcf, fullfile(folderSave, ['CoeffCount' suffix '.pdf']));