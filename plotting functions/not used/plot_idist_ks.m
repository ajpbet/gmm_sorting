function plot_idist_ks(sorted_idist, ind_idist, idist_kmatch_lSel, ks_out_full, all_ks, lenKs, folderName, channelNum)
    fig_idist = figure;
    r_sorted_idist = flip(sorted_idist);
    r_ind_idist = flip(ind_idist);
    plot(r_sorted_idist)
    title("idist vals per coeff")


    hold on
    for k = 1:length(r_sorted_idist)
        if k <= idist_kmatch_lSel  % Below x=13 - use dots
            plot(k, r_sorted_idist(k), 'bo', 'MarkerSize', 6, 'MarkerFaceColor','b');
        else  % Above or at x=13 - use x markers
            plot(k, r_sorted_idist(k), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end
    for k = 1:length(r_sorted_idist)
       text(k+0.1,r_sorted_idist(k)+0.1,num2str(r_ind_idist(k)),"Color",'b'); 
    end
    hold on;
    ks_d = flip(ks_out_full);
    ks_y = flip(all_ks);
    yyaxis right;
    plot(ks_y)
    hold on
    
    % Add scatter points with different markers based on x-position
    for k = 1:length(ks_y)
        if k <= lenKs  % Below x=13 - use dots
            plot(k, ks_y(k), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        else  % Above or at x=13 - use x markers
            plot(k, ks_y(k), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end
    
    % Add text labels (your existing code)
    for k = 1:length(ks_y)
       text(k-0.2, ks_y(k)-0.05*range(ks_y), num2str(ks_d(k)), "Color", 'r', ...
            'FontSize', 8, 'HorizontalAlignment', 'left'); 
    end
    
    % Add vertical line at x=13
    xline(lenKs, '--', 'Color', 'red', 'LineWidth', 0.5, 'Alpha', 0.7);
    xline(idist_kmatch_lSel, '--', 'Color', 'blue', 'LineWidth', 0.5, 'Alpha', 0.7);

    % Match types (colors)
    h_match(1) = plot(NaN, NaN, 'b-', 'DisplayName', 'IDIST');          % blue line
    h_match(2) = plot(NaN, NaN, 'Color', [1 0.5 0], 'LineStyle', '-', ...
                      'DisplayName', 'K');                               % orange line
    h_match(3) = scatter(NaN, NaN, 80, 'k', 'filled', 'DisplayName', 'Selected');  % filled dot
    h_match(4) = scatter(NaN, NaN, 80, 'kx', 'DisplayName', 'Not selected');       % X marker
    
    legend(h_match, 'Location', 'best');

    hold off;
    filename_IDKScoeff = fullfile(folderName,sprintf('ch%s_IDKScoeff.png', channelNum));
    exportgraphics(fig_idist, filename_IDKScoeff, 'Resolution', 300);
end