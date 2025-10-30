function idistKmatch_vKv(idist_kmatch,kDist_vec,poly_match_idistK,polyID,g,pg,xg,idist_select, ...
    idist_kmatch_lSel, channelNum, folderName)
    fig_idistMatch_kv = figure;
    title("top idist vs. kv")
    grid on
    hold on
    % For each coefficient, plot IDist (x) vs KS (y)
    for k = 1:length(idist_kmatch)    
        [pk_vals, pk_locs, pk_widths] = findpeaks(pg{k}, xg{k}, 'WidthReference', 'halfheight');
        if ~isempty(idist_kmatch(k))&&(idist_kmatch(k)~=0)
            medDist = idist_kmatch(k);
            gaussSelect = poly_match_idistK(k);
            gauss_idx = find(polyID{k} == gaussSelect);

            kv = kj_mat{k}(gauss_idx);
            % Find the IDist value for this coefficient number
            inRangeIdx = pk_locs >= gauss_low & pk_locs <= gauss_upp;
            hasPeakInRange = any(inRangeIdx);


            % Determine marker based on top 13 KS
            if ismember(k,idist_select(:,1))  % Since ks_y is sorted, first 13 are top
                % Top 13 KS - use dot
                marker = 'o';
                if hasPeakInRange
                    col = 'g';
                else
                    col = 'r';
                end
            else
                % Not top 13 KS - use x
                marker = 'x';
                if hasPeakInRange
                    col = 'g';
                else
                    col = 'r';
                end
                marker_face = 'none';
            end
            
            % Plot the point
            scatter(medDist, kv, 80, col, marker, 'LineWidth', 1.5, ...
                    'MarkerFaceColor', col);
            
            % Label with coefficient number
            text(medDist, kv+0.4, num2str(k), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
        end
    end
    xline(idist_select(idist_kmatch_lSel,2), '--', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'IDist placeholder');
    h_m(1) = scatter(NaN, NaN, 80, 'g', 'o', 'filled', 'DisplayName', 'near peak');
    h_m(2) = scatter(NaN, NaN, 80, 'r', 'o', 'filled', 'DisplayName', 'no peak');
    h_mk(1) = scatter(NaN, NaN, 80, 'k', 'o', 'filled', 'DisplayName', 'idist selected');
    h_mk(2) = scatter(NaN, NaN, 80, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'idist not selected');
    h_mk(3) = plot(NaN, NaN, '--', 'Color', 'b', 'DisplayName', 'idist knee');
   
    legend([h_m, h_mk], 'Location', 'best');

    hold off;
    xlabel('idistKmatch v kVal Values');
    ylabel('K Values'); 
    filename_idistMatch_kv = fullfile(folderName,sprintf('ch%s_idistKmatch_kv.png', channelNum));
    % exportgraphics(fig_idistMatch_kv, filename_idistMatch_kv, 'Resolution', 300);

end