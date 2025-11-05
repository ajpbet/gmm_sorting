function plotKSCoeffCount(ks_coeff_all, select_all, unitNumbers,plotTitle,graphNum)
% PLOTKSCOEFFCOUNT plots number of coefficients per unit with shading and jitter
% Inputs:
%   ks_coeff_all : cell array of ks coefficients per unit
%   select_all   : cell array of select_all coefficients per unit
%   unitNumbers  : vector of unit numbers
%   graphNum     : integer (1-5) to define suffix for saving
%   plotTitle    : string, title for the plot

numUnits = numel(unitNumbers);
figure; hold on;

% Determine max count for axis limitsks_counts  = cellfun(@numel, ks_coeff_all);
ks_counts  = cellfun(@numel, ks_coeff_all);
sel_counts = cellfun(@(x) numel(unique(x)), select_all);  % unique counts
yMax = ceil(max([ks_counts; sel_counts]) * 1.2);          % take max, scale 20%


% Background shading for alternating units
for i = 1:numUnits
    unit = unitNumbers(i);
    if mod(unit,2) == 0
        xPatch = [unit-0.5, unit+0.5, unit+0.5, unit-0.5];
        yPatch = [0, 0, 100, 100];
        patch(xPatch, yPatch, [0.95 0.95 0.95], 'EdgeColor','none'); % light gray
    end
end


% Plot points with jitte
hKS = [];
hSel = [];
jitterAmp = 0.5;  % +/- 0.1
for i = 1:numUnits
    x_jitter = unitNumbers(i) + (rand()*jitterAmp - jitterAmp/2);

    % KS coefficients
    ks_count = numel(ks_coeff_all{i});
    if isempty(hKS)
        hKS = plot(x_jitter, ks_count, 'o', 'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerEdgeColor','k');
    else
        plot(x_jitter, ks_count, 'o', 'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerEdgeColor','k');
    end

    % select_all coefficients
    sel_count = numel(unique(select_all{i}));
    if isempty(hSel)
        hSel = plot(x_jitter, sel_count, 's', 'MarkerFaceColor', [0.8 0.4 0.2], 'MarkerEdgeColor','k');
    else
        plot(x_jitter, sel_count, 's', 'MarkerFaceColor', [0.8 0.4 0.2], 'MarkerEdgeColor','k');
    end
end

% Formatting 
xlabel('Unit Number');
ylabel('Number of Coefficients');
title(plotTitle);
xlim([min(unitNumbers)-0.5, max(unitNumbers)+0.5]);
ylim([0 yMax]);
grid on;
legend([hKS, hSel], {'KS Coefficients','GMM Coefficients'}, 'Location','best');

% --- Save figure ---
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

saveas(gcf, fullfile(folderSave, ['CoeffCount' suffix '.png']));
end
