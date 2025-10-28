  %% best idist and its respective k val
    fig_idist_kv = figure;
    title("idist v kv")
    grid on
    hold on
    % For each coefficient, plot IDist (x) vs KS (y)
    for k = 1:length(r_sorted_idist)
        coeff_num = r_ind_idist(k);  % Coefficient number from KS ordering
        gaussWinner = poly_match_idistK(coeff_num);
        newGaussIdx = find(polyID{coeff_num} == gaussWinner);


        kval = kj_mat{coeff_num}(newGaussIdx);
        if isempty(newGaussIdx)
            kval = 0;
        end
        % Find the IDist value for this coefficient number
        idist_val = idist_kmatch(coeff_num);
        
        % Determine marker based on top 13 KS
        if k <= idist_kmatch_lSel  % Since ks_y is sorted, first 13 are top
            % Top 13 KS - use dot
            marker = 'o';
            marker_face = color_idistk{coeff_num};
        else
            % Not top 13 KS - use x
            marker = 'x';
            marker_face = 'none';
        end
        
        % Plot the point
        scatter(idist_val, kval, 80, color_idistk{coeff_num}, marker, 'LineWidth', 1.5, ...
                'MarkerFaceColor', marker_face);
        
        % Label with coefficient number
        text(idist_val, kval+0.01, num2str(coeff_num), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
    end

    h_k(1) = scatter(NaN, NaN, 80, 'k', 'o', 'filled', 'DisplayName', 'idist selected');
    h_k(2) = scatter(NaN, NaN, 80, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'idist not selected');

    legend(h_k, 'Location', 'best');
    hold off;
    xlabel('IDist Values');
    ylabel('K Values'); 
    filename_idist_kv = fullfile(folderName,sprintf('ch%s_id_v_kv.png', channelNum));
  %  exportgraphics(fig_idist_kv, filename_idist_kv, 'Resolution', 300);