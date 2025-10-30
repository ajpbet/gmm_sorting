function plotMedDistVsKv(medDist_vec, kDist_vec, medDist_select, kDist_select, select_gauss, folderName, channelNum, dataSetFlag)
% PLOTMEDDISTVSKV: Plot median distance vs. k-value with threshold lines and selected Gaussians
%
% Inputs:
%   medDist_vec    - [medDist, coeff#, gauss#]
%   kDist_vec      - [k_value, coeff#, gauss#]
%   medDist_select - median distance intersection points
%   kDist_select   - k-value intersection points
%   select_gauss   - [medDist, k_value, coeff#, gauss#] of selected Gaussians
%   folderName     - folder path to save figure
%   channelNum     - channel identifier (numeric or string)
%   dataSetFlag    - 1 = "all peaks", 2 = "pks excluded"

    % Determine title and file suffix based on dataSetFlag
    switch dataSetFlag
        case 1
            figTitle = 'Median Distance vs K-value (All Peaks)';
            fileSuffix = 'all_peaks';
        case 2
            figTitle = 'Median Distance vs K-value (Pks Excluded)';
            fileSuffix = 'pks_excluded';
        otherwise
            error('dataSetFlag must be 1 or 2');
    end

    figure('Visible','on'); hold on; grid on; box on;

    % --- Compute Q2/Q4 lines based on intersection
    x_inter = medDist_select(end,1);
    y_inter = kDist_select(end,1);

    % Limits
    maxKv = max(kDist_vec(:,1));
    minKv = min(kDist_vec(:,1));
    maxMedDist = max(medDist_vec(:,1));
    minMedDist = min(medDist_vec(:,1));

    % Q2 line (intersect → minMedDist, maxKv)
    x_Q2 = [x_inter, minMedDist];
    y_Q2 = [y_inter, maxKv];
    plot(x_Q2, y_Q2, 'k--', 'LineWidth', 1.2);

    % Q4 line (intersect → maxMedDist, minKv)
    x_Q4 = [x_inter, maxMedDist];
    y_Q4 = [y_inter, minKv];
    plot(x_Q4, y_Q4, 'k--', 'LineWidth', 1.2);

    % --- Plot all points with color/marker scheme
    colorMap = struct('g',[0 1 0],'b',[0 0 1],'r',[1 0 0]);
    
    for k = 1:size(medDist_vec,1)
        medDist = medDist_vec(k,1);
        kv      = kDist_vec(kDist_vec(:,2)==medDist_vec(k,2) & kDist_vec(:,3)==medDist_vec(k,3),1);
        coeff_num = medDist_vec(k,2);
        gauss_idx = medDist_vec(k,3);
    
        % Determine color/marker
        if ismember([coeff_num, gauss_idx], select_gauss(:,3:4),'rows')
            col = 'g'; marker = 'o'; % selected Gaussian
        elseif ismember(coeff_num, select_gauss(:,3))
            col = 'b'; marker = 'o'; % same coefficient but not selected
        else
            col = 'r'; marker = 'x'; % unrelated
        end
    
        if marker == 'o'
            scatter(medDist, kv, 80, colorMap.(col), marker, 'LineWidth', 1.5, 'MarkerFaceColor', colorMap.(col));
        else
            scatter(medDist, kv, 80, colorMap.(col), marker, 'LineWidth', 1.5); % x marker unfilled
        end
    
        text(medDist, kv+0.4, num2str(coeff_num), 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', 'FontSize',8, 'Color','k', 'FontWeight','bold');
    end

    % Intersection point
    plot(x_inter, y_inter, 'ro', 'MarkerSize',8, 'LineWidth',1.5);

    % --- Plot reference xline and yline (from intersection selection)
    xline_val = medDist_select(end,1); % last value in medDist_select
    yline_val = kDist_select(end,1);   % last value in kDist_select
    
    h_xline = xline(xline_val, '--', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Reference medDist');
    h_yline = yline(yline_val, '--', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Reference kValue');
    
    % Add to legend (after scatter points)

    
    % --- Legend
    h(1) = scatter(NaN, NaN, 80, 'g', 'o', 'filled', 'DisplayName','Selected Gaussian');
    h(2) = scatter(NaN, NaN, 80, 'b', 'o', 'filled', 'DisplayName','Same Coefficient');
    h(3) = scatter(NaN, NaN, 80, 'r', 'x', 'LineWidth',1.5, 'DisplayName','Unrelated');
    h(4) = plot(NaN, NaN, 'k--', 'LineWidth',1.2, 'DisplayName','Threshold lines');

    legend([h(1:4), h_xline, h_yline],'Location','best');

    xlabel('Median Distance');
    ylabel('K-value');
    title(figTitle);

    % --- Save figure
    if isnumeric(channelNum)
        channelNum = num2str(channelNum);
    end
    if ~exist(folderName,'dir')
        mkdir(folderName);
    end
    filename = fullfile(folderName, sprintf('ch%s_medD_v_kv_%s.png', channelNum, fileSuffix));
    exportgraphics(gcf, filename, 'Resolution',300);

    hold off;
end
