function plotKSCoeffMismatch(ks_coeff_all, select_all, unitNumbers, plotTitle, graphNum)
% PLOTKSCOEFFMISMATCH plots percent mismatch between ks_coeff and select_all gaussians
%   - % of ks_coeff not in select_all
%   - % of select_all not in ks
%   Also adds a paired sign test to the scatter plot and generates a separate
%   box plot of the mismatch distributions.
%   Saves figures to results/mismatches with appropriate suffix.
%
% Inputs:
%   ks_coeff_all : cell array or matrix of reference coefficients (per unit)
%   select_all   : cell array or matrix of selected coefficients (per unit)
%   unitNumbers  : vector of unit numbers (e.g., 1,2,...,N)
%   plotTitle    : string, title for the plot
%   graphNum     : integer 1–5, used for suffix (_0pct, _1pct, etc.)

% --- Original Scatter Plot ---
figure; hold on;

% Initialize arrays to store all percentage values for stats
all_pct_ks = [];
all_pct_sel = [];

numUnits = numel(unitNumbers);
for i = 1:numUnits
    unit = unitNumbers(i);
    if mod(unit,2) == 0
        xPatch = [unit-0.5, unit+0.5, unit+0.5, unit-0.5];
        yPatch = [0, 0, 100, 100];
        patch(xPatch, yPatch, [0.95 0.95 0.95], 'EdgeColor','none'); % light gray
    end
end

% Initialize handle arrays
hKS  = plot(nan, nan, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], 'MarkerEdgeColor','k');
hSel = plot(nan, nan, 's', 'MarkerFaceColor',[0.8 0.4 0.2], 'MarkerEdgeColor','k');

for i = 1:numUnits
    ks = ks_coeff_all{i};
    selCell = unique(select_all{i});
    
    if ~iscell(selCell)
        selCell = {selCell};
    end
    
    for t = 1:numel(selCell)
        sel = selCell{t};
        ksVec = ks(:);
        selVec = sel(:);
        
        % Use rounding to account for potential precision issues
        ksVecRound = round(ksVec, 6);
        selVecRound = round(selVec, 6);
        overlap = intersect(ksVecRound, selVecRound);
        
        % Calculate total as the size of the union
        total = numel(ksVecRound) + numel(selVecRound) - numel(overlap);
        
        if isempty(ksVec) || total == 0
            pct_ks = 0;
        else
            pct_ks = (numel(ksVecRound) - numel(overlap)) / total * 100;
        end
        
        if isempty(selVec) || total == 0
            pct_sel = 0;
        else
            pct_sel = (numel(selVecRound) - numel(overlap)) / total * 100;
        end
        
        % Store values for stats
        all_pct_ks = [all_pct_ks; pct_ks];
        all_pct_sel = [all_pct_sel; pct_sel];
        
        x_jitter = unitNumbers(i) + (rand()*0.5 - 0.1);
        
        plot(x_jitter, pct_ks, 'o', 'MarkerFaceColor',[0.2 0.6 0.8], 'MarkerEdgeColor', 'k');
        plot(x_jitter, pct_sel, 's', 'MarkerFaceColor', [0.8 0.4 0.2], 'MarkerEdgeColor', 'k');
    end
end

xlabel('Unit Number');
ylabel('Percent Mismatch (%)');
ylim([0 100]);
xlim([min(unitNumbers)-0.5, max(unitNumbers)+0.5]);
title(plotTitle);
grid on;
box on;
legend([hKS, hSel], {'% KS not in GMM','% GMM not in KS'}, 'Location','best');

% --- Add Paired Sign Test to Scatter Plot ---
if ~isempty(all_pct_ks)
    % Run the paired sign test
    [p, ~, ~] = signtest(all_pct_ks, all_pct_sel);
    textStr = sprintf('Paired Sign Test: p = %.3e', p);
else
    textStr = 'Paired Sign Test: No data';
end
% Add text to the top-left corner
ax = gca;
text(ax.XLim(1) + 0.02*range(ax.XLim), ax.YLim(2) - 0.02*range(ax.YLim), ...
     textStr, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
     'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k');


% --- Save the Scatter Plot ---
switch graphNum
    case 1, suffix = '_0pct';
    case 2, suffix = '_1pct';
    case 3, suffix = '_2p5pct';
    case 4, suffix = '_5pct';
    case 5, suffix = '_10pct';
    otherwise, suffix = '_spikeMatch';
end
saveFolder = fullfile('results', 'mismatches');
if ~exist(saveFolder,'dir'), mkdir(saveFolder); end
fileNamePNG = ['ks_mismatch' suffix '.pdf']; % Original code used .pdf
saveas(gcf, fullfile(saveFolder, fileNamePNG));

% --- Create and Save Additional Box Plot ---
if ~isempty(all_pct_ks)
    figure; % Create a new figure
    
    data_for_boxplot = [all_pct_ks, all_pct_sel];
    labels_for_boxplot = {'% KS not in GMM', '% GMM not in KS'};
    
    boxplot(data_for_boxplot, 'Labels', labels_for_boxplot);
    
    ylabel('Percent Mismatch (%)');
    % Use a multi-line title for clarity
    title({plotTitle; 'Distribution of Mismatch Percentages (All Units & Tests)'});
    grid on;
    box on;
    
    % Save the box plot to a separate file
    fileNameBoxPlotPDF = ['ks_mismatch' suffix '_boxplot.pdf'];
    saveas(gcf, fullfile(saveFolder, fileNameBoxPlotPDF));
end

end % End of function