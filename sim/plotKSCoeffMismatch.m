function plotKSCoeffMismatch(ks_coeff_all, select_all, unitNumbers, plotTitle, graphNum)
% PLOTKSCOEFFMISMATCH plots percent mismatch between ks_coeff and select_all gaussians
%   - % of ks_coeff not in select_all
%   - % of select_all not in ks
%   Saves figure to results/mismatches with appropriate suffix.
%
% Inputs:
%   ks_coeff_all : cell array or matrix of reference coefficients (per unit)
%   select_all   : cell array or matrix of selected coefficients (per unit)
%   unitNumbers  : vector of unit numbers (e.g., 1,2,...,N)
%   plotTitle    : string, title for the plot
%   graphNum     : integer 1–5, used for suffix (_0pct, _1pct, etc.)

figure; hold on;
numUnits = numel(unitNumbers);

% --- Background shading for alternating units ---
for i = 1:numUnits
    unit = unitNumbers(i);
    if mod(unit,2) == 0
        xPatch = [unit-0.5, unit+0.5, unit+0.5, unit-0.5];
        yPatch = [0, 0, 100, 100];
        patch(xPatch, yPatch, [0.95 0.95 0.95], 'EdgeColor','none'); % light gray
    end
end

% Initialize handle arrays
hKS = [];
hSel = [];

% --- Plot points with small random jitter ---
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
        
        if isempty(ksVec)
            pct_ks = 0;
        else
            overlap = intersect(round(ksVec,6), round(selVec,6));
            pct_ks = (numel(ksVec)-numel(overlap))/numel(ksVec)*100;
        end
        
        if isempty(selVec)
            pct_sel = 0;
        else
            overlap = intersect(round(ksVec,6), round(selVec,6));
            pct_sel = (numel(selVec)-numel(overlap))/numel(selVec)*100;
        end
        
        x_jitter = unitNumbers(i) + (rand()*0.5 - 0.1);
        
        % Store handles on first plot only for legend
        if isempty(hKS)
            hKS = plot(x_jitter, pct_ks, 'o', 'MarkerFaceColor', [0.85 0.3 0.3], 'MarkerEdgeColor', 'k');
        else
            plot(x_jitter, pct_ks, 'o', 'MarkerFaceColor', [0.85 0.3 0.3], 'MarkerEdgeColor', 'k');
        end
        
        if isempty(hSel)
            hSel = plot(x_jitter, pct_sel, 's', 'MarkerFaceColor', [0.3 0.6 0.3], 'MarkerEdgeColor', 'k');
        else
            plot(x_jitter, pct_sel, 's', 'MarkerFaceColor', [0.3 0.6 0.3], 'MarkerEdgeColor', 'k');
        end
    end
end
xlabel('Unit Number');
ylabel('Percent Mismatch (%)');
ylim([0 100]);
xlim([min(unitNumbers)-0.5, max(unitNumbers)+0.5]);
title(plotTitle);
grid on;
box on;
legend([hKS, hSel], {'% KS not in select\_all','% GMM not in KS'}, 'Location','best');

% --- Determine suffix based on graphNum ---
switch graphNum
    case 1, suffix = '_0pct';
    case 2, suffix = '_1pct';
    case 3, suffix = '_2p5pct';
    case 4, suffix = '_5pct';
    case 5, suffix = '_10pct';
    otherwise, suffix = sprintf('_graph%d', graphNum);
end

% --- Save figure ---
saveFolder = fullfile('results', 'mismatches');
if ~exist(saveFolder,'dir'), mkdir(saveFolder); end

fileNamePNG = ['ks_mismatch' suffix '.png'];

saveas(gcf, fullfile(saveFolder, fileNamePNG));

end
