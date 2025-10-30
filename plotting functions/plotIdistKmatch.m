function plotIdistKmatch(all_ks,ks_out_full, lenKs, idist_kmatch, color_idistk, ks_out, idist_select, folderName, channelNum)
% PLOTIDISTKMATCH: Plot IDist vs KS values highlighting top KS coefficients
%
% Inputs:
%   ks_y          - Sorted KS values
%   ks_d          - Coefficient numbers corresponding to sorted KS values
%   lenKs         - Number of top KS coefficients to highlight (e.g., 13)
%   idist_kmatch  - Vector mapping coefficient number to IDist value
%   color_idistk  - Cell array of colors per coefficient
%   ks_out        - Matrix of KS results (for cutoff line)
%   idist_select  - Matrix of selected IDist values (for cutoff line)
%   folderName    - Output folder path for saving the figure
%   channelNum    - Channel identifier (string or numeric)
%
% Example:
%   plotIdistKmatch(ks_y, ks_d, 13, idist_kmatch, color_idistk, ks_out, idist_select, 'output', '2');

    ks_d = flip(ks_out_full);
    ks_y = flip(all_ks);
    % Create figure
    fig1st_idistKmatch = figure('Visible', 'on'); % Change to 'off' for silent plotting
    hold on;

    % Plot each coefficient point
    for k = 1:length(ks_y)
        coeff_num = ks_d(k);  % Coefficient index in KS order

        % Defensive: ensure coefficient index is valid
        if coeff_num > numel(idist_kmatch)
            continue;
        end

        % Find IDist value for this coefficient
        idist_val = idist_kmatch(coeff_num);

        % Determine marker and fill style
        if k <= lenKs
            marker = 'o';  % Top KS - filled circle
            marker_face = color_idistk{coeff_num};
        else
            marker = 'x';  % Others - cross
            marker_face = 'none';
        end

        % Plot point
        scatter(idist_val, ks_y(k), 80, color_idistk{coeff_num}, marker, ...
                'LineWidth', 1.5, 'MarkerFaceColor', marker_face);

        % Label coefficient number
        text(idist_val, ks_y(k) + 0.01, num2str(coeff_num), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
    end

    % Axes and title
    xlabel('IDist Values');
    ylabel('KS Values');
    title('IDist vs KS Values - highest idist with top 3 k');
    grid on;

    % Add threshold lines
    if ~isempty(ks_out)
        yline(min(ks_out(:,2)), '--', 'Color', [1 0.5 0], ...
            'LineWidth', 1.5, 'DisplayName', 'KS cutoff');
    end

    if ~isempty(idist_select)
        xline(min(idist_select(:,2)), '--', 'Color', 'b', ...
            'LineWidth', 1.5, 'DisplayName', 'IDist cutoff');
    end

    % Create legend handles
    h_match(1) = scatter(NaN, NaN, 80, 'm', 'o', 'filled', 'DisplayName', '1st mDist');
    h_match(2) = scatter(NaN, NaN, 80, 'g', 'o', 'filled', 'DisplayName', '2nd mDist');
    h_match(3) = scatter(NaN, NaN, 80, 'b', 'o', 'filled', 'DisplayName', '3rd mDist');
    h_match(4) = scatter(NaN, NaN, 80, 'r', 'o', 'filled', 'DisplayName', '>3 mDist');
    h_match(5) = plot(NaN, NaN, '--', 'Color', [1 0.5 0], 'DisplayName', 'KS cutoff');
    h_match(6) = plot(NaN, NaN, '--', 'Color', 'b', 'DisplayName', 'IDist cutoff');

    h_ks(1) = scatter(NaN, NaN, 80, 'k', 'o', 'filled', 'DisplayName', 'KS selected');
    h_ks(2) = scatter(NaN, NaN, 80, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'KS not selected');

    legend([h_match, h_ks], 'Location', 'best');
    hold off;

    % Save figure
    if isnumeric(channelNum)
        channelNum = num2str(channelNum);
    end
    filename_idistKmatch = fullfile(folderName, sprintf('ch%s_idistKmatch.png', channelNum));
    exportgraphics(fig1st_idistKmatch, filename_idistKmatch, 'Resolution', 300);
end
