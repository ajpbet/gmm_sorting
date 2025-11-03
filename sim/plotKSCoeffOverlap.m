function plotKSCoeffOverlap(ks_coeff_all, select_all, unitNumbers, plotTitle, graphNum)
% PLOTKSCOEFFOVERLAP plots percent overlap between ks_coeff and select_all gaussians
% and saves the figure with a unique identifier based on percentage.

figure; hold on;

numUnits = numel(unitNumbers);

% --- Background shading for alternating units ---
for i = 1:numUnits
    unit = unitNumbers(i);
    if mod(unit,2) == 0 % shade even units
        xPatch = [unit-0.5, unit+0.5, unit+0.5, unit-0.5];
        yPatch = [0, 0, 100, 100];
        patch(xPatch, yPatch, [0.95 0.95 0.95], 'EdgeColor','none'); % light gray
    end
end

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
        
        if isempty(ksVec) || isempty(selVec)
            pct = 0;
        else
            overlap = intersect(round(ksVec,6), round(selVec,6));
            pct = numel(overlap)/numel(ksVec)*100;
        end
        
        % Random jitter in x-direction (±0.1)
        x_jitter = unitNumbers(i) + (rand()*0.5 - 0.1);
        plot(x_jitter, pct, 'o', 'MarkerFaceColor', [0.2 0.6 0.8], 'MarkerEdgeColor', 'k');
    end
end

xlabel('Unit Number');
ylabel('Percent Overlap (%)');
ylim([0 100]);
xlim([min(unitNumbers)-0.5, max(unitNumbers)+0.5]);
title(plotTitle);
grid on;
box on;

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
saveFolder = fullfile('results', 'matches');
if ~exist(saveFolder,'dir'), mkdir(saveFolder); end

%fileNameFig = ['ks_overlap' suffix '.fig'];
fileNamePNG = ['ks_overlap' suffix '.png'];

%savefig(fullfile(saveFolder, fileNameFig));
saveas(gcf, fullfile(saveFolder, fileNamePNG));

end
