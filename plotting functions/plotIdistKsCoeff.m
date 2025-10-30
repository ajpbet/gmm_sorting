function plotIdistKsCoeff(sorted_idist, ind_idist, idist_kmatch_lSel, ks_out_full, all_ks, lenKs, folderName, channelNum)
% PLOTIDISTKSCOEFF: Plot sorted IDist and KS values across coefficients
%
% Inputs:
%   sorted_idist        - Vector of IDist values (ascending or descending order)
%   ind_idist           - Corresponding coefficient indices for IDist values
%   idist_kmatch_lSel   - X-position (index) cutoff for IDist selection
%   ks_out_full         - Vector or matrix containing KS coefficient indices
%   all_ks              - Vector of KS values (same length as ks_out_full)
%   lenKs               - X-position (index) cutoff for KS selection
%   folderName          - Output folder path for saving the figure
%   channelNum          - Channel identifier (string or numeric)
%
% Outputs:
%   filename_IDKScoeff  - Full path to the saved PNG figure
%
% Example:
%   plotIdistKsCoeff(sorted_idist, ind_idist, 13, ks_out_full, all_ks, 13, 'output', '5');

    % Reverse (flip) order for plotting
    r_sorted_idist = flip(sorted_idist);
    r_ind_idist = flip(ind_idist);



    % Create figure
    fig_idist = figure('Visible', 'on'); % Change to 'off' for batch mode
    hold on
    grid on
    title('IDist and KS Values per Coefficient');
    xlabel('Coefficient Index (Sorted)');
    ylabel('IDist Value');
    plot(r_sorted_idist, 'b-', 'DisplayName', 'IDIST');
    
 
    % Plot IDist points and markers
    for k = 1:length(r_sorted_idist)
        if k <= idist_kmatch_lSel
            plot(k, r_sorted_idist(k), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
        else
            plot(k, r_sorted_idist(k), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end

    % Label coefficient numbers for IDist
    for k = 1:length(r_sorted_idist)
        text(k + 0.1, r_sorted_idist(k) + 0.1, num2str(r_ind_idist(k)), 'Color', 'b', ...
             'FontSize', 8, 'HorizontalAlignment', 'left');
    end

    % Flip KS arrays for alignment
    ks_d = flip(ks_out_full);
    ks_y = flip(all_ks);

    % Plot KS values on right axis
    yyaxis right;
    plot(ks_y, 'Color', [1 0.5 0], 'LineStyle', '-', 'DisplayName', 'K');
    ylabel('KS Value');

    % Plot KS points with marker distinction
    for k = 1:length(ks_y)
        if k <= lenKs
            plot(k, ks_y(k), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        else
            plot(k, ks_y(k), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end

    % Label coefficient numbers for KS
    for k = 1:length(ks_y)
        text(k - 0.2, ks_y(k) - 0.05 * range(ks_y), num2str(ks_d(k)), 'Color', 'r', ...
             'FontSize', 8, 'HorizontalAlignment', 'left');
    end

    % Add reference lines
    xline(lenKs, '--', 'Color', 'red', 'LineWidth', 0.5, 'Alpha', 0.7, 'DisplayName', 'KS cutoff');
    xline(idist_kmatch_lSel, '--', 'Color', 'blue', 'LineWidth', 0.5, 'Alpha', 0.7, 'DisplayName', 'IDist cutoff');

    % Build legend
    h_match(1) = plot(NaN, NaN, 'b-', 'DisplayName', 'IDIST');
    h_match(2) = plot(NaN, NaN, 'Color', [1 0.5 0], 'LineStyle', '-', 'DisplayName', 'K');
    h_match(3) = scatter(NaN, NaN, 80, 'k', 'filled', 'DisplayName', 'Selected');
    h_match(4) = scatter(NaN, NaN, 80, 'kx', 'DisplayName', 'Not selected');

    legend(h_match, 'Location', 'best');
    hold off;

    % Format output filename and save
    if isnumeric(channelNum)
        channelNum = num2str(channelNum);
    end
    filename_IDKScoeff = fullfile(folderName, sprintf('ch%s_IDKScoeff.png', channelNum));
    exportgraphics(fig_idist, filename_IDKScoeff, 'Resolution', 300);

    % Return output path
    if nargout > 0
        varargout{1} = filename_IDKScoeff;
    end
end
