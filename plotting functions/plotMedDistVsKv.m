function plotMedDistVsKv(pg, xg, medDistVector, polyID, g, kDist_vec, medD_sel_noPk, medD_lSel_noPk, ...
     k_sel_NoPk,k_lSel_NoPk,folderName, channelNum)
% PLOTMEDDISTVSKV: Plot median distance vs. k-value, highlighting top IDist selections and peak ranges
%
% Inputs:
%   pg, xg                - Cell arrays of Gaussian PDFs and corresponding x-values
%   medDist_sort          - Cell array of median distance values per coefficient
%   medDist_sortIdx       - Cell array of sorted index selections per coefficient
%   polyID                - Cell array mapping Gaussian index to polynomial ID
%   g                     - Cell array of Gaussian mixture objects (with fields mu and Sigma)
%   kj_mat                - Cell array of K-values per Gaussian component
%   idist_select          - Matrix of selected IDist indices and distances (cols: index, distance)
%   idist_kmatch_lSel     - Indices used to mark reference x-lines
%   folderName            - Output folder path for saving the figure
%   channelNum            - Channel identifier (string or numeric)
%
% Example:
%   plotMedDistVsKv(pg, xg, medDist_sort, medDist_sortIdx, polyID, g, kj_mat, idist_select, idist_kmatch_lSel, 'output', '5');

    fig_meddist_kv = figure('Visible', 'on'); % change to 'off' if you want silent plotting
    title("meddist vs. kv")
    grid on
    hold on

    % Loop over all coefficients
    for k = 1:length(medDistVector(:,1))
        medDist = medDistVector(k,1);
        coeff_num = medDistVector(k,2);
        gauss_idx = medDistVector(k,3);
        
        kv_coeffs = kDist_vec(kDist_vec(:,2)==coeff_num,:);
        kv = kv_coeffs(kv_coeffs(:,3)==gauss_idx); % Retrieve the k-value for the current coefficient

        % Determine marker and color based on IDist selection and peak
        if ismember(coeff_num, medD_sel_noPk(:,2))
            marker = 'o';
            col = 'g';
        elseif ismember(coeff_num,k_sel_NoPk(:,2))
            marker = 'o';
            col = 'b';
        else
            marker = 'x';
            col = 'r';
            
        end

        % Plot scatter point
        scatter(medDist, kv, 80, col, marker, 'LineWidth', 1.5, 'MarkerFaceColor', col);

        % Label coefficient number
        text(medDist, kv + 0.4, num2str(coeff_num), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
                'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
        
    end

    % Add IDist reference line
    xline(medD_sel_noPk(medD_lSel_noPk,1), '--', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'IDist placeholder');
    yline(k_sel_NoPk( ...
        k_lSel_NoPk,1), '--', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'IDist placeholder')
    
    % Create legend handles
    h_m(1) = scatter(NaN, NaN, 80, 'g', 'o', 'filled', 'DisplayName', 'near peak');
    h_m(2) = scatter(NaN, NaN, 80, 'r', 'o', 'filled', 'DisplayName', 'no peak');
    h_mk(1) = scatter(NaN, NaN, 80, 'k', 'o', 'filled', 'DisplayName', 'idist selected');
    h_mk(2) = scatter(NaN, NaN, 80, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'idist not selected');
    h_mk(3) = plot(NaN, NaN, '--', 'Color', 'b', 'DisplayName', 'idist knee');

    legend([h_m, h_mk], 'Location', 'best');
    hold off;

    xlabel('MedDist Values');
    ylabel('K Values');

    % Save figure
    if isnumeric(channelNum)
        channelNum = num2str(channelNum);
    end
    filename_meddist_kv = fullfile(folderName, sprintf('ch%s_medD_v_kv.png', channelNum));
    %exportgraphics(fig_meddist_kv, filename_meddist_kv, 'Resolution', 300);
end
