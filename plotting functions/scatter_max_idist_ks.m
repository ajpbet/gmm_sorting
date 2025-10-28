function scatter_max_idist_ks(ks_d,ks_y,lenKs)
    fig_maxIDvKS = figure;
    hold on;
    
    % Top 13 KS coefficients are the first 13 since ks_y is sorted
    ks_top_coeffs = ks_d(1:lenKs);
    
    % For each coefficient, plot IDist (x) vs KS (y)
    for k = 1:length(ks_y)
        coeff_num = ks_d(k);  % Coefficient number from KS ordering
        
        % Find the IDist value for this coefficient number
        idist_val = idistVals(coeff_num);
        
        % Determine marker based on top 13 KS
        if k <= lenKs  % Since ks_y is sorted, first 13 are top
            % Top 13 KS - use dot
            marker = 'o';
            marker_face = colors{coeff_num};
        else
            % Not top 13 KS - use x
            marker = 'x';
            marker_face = 'none';
        end
        
            % Plot the point
        scatter(idist_val, ks_y(k), 80, colors{coeff_num}, marker, 'LineWidth', 1.5, ...
                'MarkerFaceColor', marker_face);
        
        % Label with coefficient number
        text(idist_val, ks_y(k)+0.01, num2str(coeff_num), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
    end
    
    xlabel('IDist Values');
    ylabel('KS Values'); 
    title('IDist vs KS Values (dots = KS select, X = Others)');
    grid on;
    
    hold on;

    % Match types (colors)
    h_match(1) = scatter(NaN, NaN, 80, 'r', 'o', 'filled', 'DisplayName', 'No matches');
    h_match(2) = scatter(NaN, NaN, 80, 'y', 'o', 'filled', 'DisplayName', '1 match');
    h_match(3) = scatter(NaN, NaN, 80, 'b', 'o', 'filled', 'DisplayName', '2 matches');
    h_match(4) = scatter(NaN, NaN, 80, 'g', 'o', 'filled', 'DisplayName', '3 matches');
    h_match(5) = scatter(NaN, NaN, 80, 'm', 'o', 'filled', 'DisplayName', 'Perfect match');
    
    % KS selection (markers)
    h_ks(1) = scatter(NaN, NaN, 80, 'k', 'o', 'filled', 'DisplayName', 'KS selected');
    h_ks(2) = scatter(NaN, NaN, 80, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'KS not selected');
    
    legend([h_match, h_ks], 'Location', 'best');
    hold off;
    filename_MaxIDvKS = fullfile(folderName,sprintf('ch%s_MaxIDvKS.png', channelNum));
    exportgraphics(fig_maxIDvKS, filename_MaxIDvKS, 'Resolution', 300);