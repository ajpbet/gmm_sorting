function plotMaxIDvKS(idistVals, all_ks, ks_out_full, colors, lenKs, folderName, channelNum)
% PLOTMAXIDVKS: Create and save a scatter plot comparing IDist and KS values
%
% INPUTS:
%   idistVals   - Vector of IDist values indexed by coefficient number
%   ks_y        - Sorted KS values (descending or ascending depending on prior step)
%   ks_d        - Coefficient indices sorted according to KS ranking
%   colors      - Cell array of RGB or color character vectors per coefficient
%   lenKs       - Number of top KS coefficients to mark (e.g., 13)
%   folderName  - Folder path to save figure
%   channelNum  - Channel number (used in filename and figure title)
%
% OUTPUT:
%   Saves figure to [folderName]/ch[channelNum]_MaxIDvKS.png
    ks_d = flip(ks_out_full);
    ks_y = flip(all_ks);
    % Create figure
    fig_maxIDvKS = figure; % Hide figure for batch processing
    hold on;
    
    % Top KS coefficients are the first lenKs since ks_y is sorted
    ks_top_coeffs = ks_d(1:lenKs);
    
    % For each coefficient, plot IDist (x) vs KS (y)
    for k = 1:length(ks_y)
        coeff_num = ks_d(k);  % Coefficient number from KS ordering
        
        % Find corresponding IDist value
        idist_val = idistVals(coeff_num);
        
        % Determine marker style for top KS
        if k <= lenKs
            marker = 'o';                      % Top KS coefficients
            marker_face = colors{coeff_num};
        else
            marker = 'x';                      % Non-top KS
            marker_face = 'none';
        end
        
        % Plot the point
        scatter(idist_val, ks_y(k), 80, colors{coeff_num}, marker, ...
            'LineWidth', 1.5, 'MarkerFaceColor', marker_face);
        
        % Label with coefficient number
        text(idist_val, ks_y(k)+0.01, num2str(coeff_num), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
    end
    
    % Axes labels and title
    xlabel('IDist Values');
    ylabel('KS Values');
    title(sprintf('IDist vs KS Values (dots = KS select, X = Others) — Ch%s', channelNum));
    grid on;
    
    % Dummy handles for legend (color → match type)
    h_match(1) = scatter(NaN, NaN, 80, 'r', 'o', 'filled', 'DisplayName', 'No matches');
    h_match(2) = scatter(NaN, NaN, 80, 'y', 'o', 'filled', 'DisplayName', '1 match');
    h_match(3) = scatter(NaN, NaN, 80, 'b', 'o', 'filled', 'DisplayName', '2 matches');
    h_match(4) = scatter(NaN, NaN, 80, 'g', 'o', 'filled', 'DisplayName', '3 matches');
    h_match(5) = scatter(NaN, NaN, 80, 'm', 'o', 'filled', 'DisplayName', 'Perfect match');
    
    % Dummy handles for KS selection markers
    h_ks(1) = scatter(NaN, NaN, 80, 'k', 'o', 'filled', 'DisplayName', 'KS selected');
    h_ks(2) = scatter(NaN, NaN, 80, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'KS not selected');
    
    legend([h_match, h_ks], 'Location', 'best');
    hold off;

    % Save figure
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    filename_MaxIDvKS = fullfile(folderName, sprintf('ch%s_MaxIDvKS.png', channelNum));
    exportgraphics(fig_maxIDvKS, filename_MaxIDvKS, 'Resolution', 300);
    

end
