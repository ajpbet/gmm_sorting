function plotSelectedVsAll(all_x, all_y, all_gauss_ids, x_inter, y_inter, trim_struct_10, selected_gauss_10pct, select_spike_match, folderName, channelNumStr, plot_and_save,ks_coeff)
% PLOTSELECTEDVSALL (Figure 3 - Modified)
% Plots LineExclusion results combined with spikeMatch redundancy filtering.
%
% Creates three categories:
% 1. Not Selected (Black 'x'): Points outside the 10% trim boundary.
% 2. Redundant (Green 'x'): Points INSIDE 10% trim, but REMOVED by spikeMatch.
% 3. Selected & Kept (Green 'o'): Points INSIDE 10% trim AND KEPT by spikeMatch.
% 4. KS Coeffs (Blue 'x' or 'o'): Points from ks_coeff list, plotted on top.
%
% MODIFIED: Added summary text box with before/after counts.
% MODIFIED: Added ks_coeff input and plotting (preserves shape).
    if ~plot_and_save
        return;
    end
    % Get ID lists for comparison
    % Set A: All points selected by 10% trim (CoeffID, GaussID)
    Set_A_IDs = selected_gauss_10pct(:, 3:4); 
    % Set B: Points KEPT by spikeMatch (subset of A)
    Set_B_IDs = select_spike_match(:, 3:4);
    
    % Initialize containers for plot points
    pts_not_selected = []; % Black 'x'
    pts_kept = [];         % Green 'o'
    pts_redundant = [];    % Green 'x'
    
    % --- NEW: KS-specific containers to preserve shape ---
    ks_pts_not_selected = []; % Blue 'x'
    ks_pts_kept = [];         % Blue 'o'
    ks_pts_redundant = [];    % Blue 'x'
    
    % Classify every single point from the original set
    for i = 1:length(all_x)
        current_id = all_gauss_ids(i, :); % Assuming [CoeffID, GaussID]
        current_pt = [all_x(i), all_y(i)];
        
        % Check if this point's coefficient is in the ks_coeff list
        is_ks = ismember(current_id(1), ks_coeff); % Assuming 1st col is CoeffID
        
        % Check if this point was inside the 10% trim boundary
        is_in_A = ismember(current_id, Set_A_IDs, 'rows');
        
        if ~is_in_A
            % 1. Not selected by lineExclusion
            if is_ks
                ks_pts_not_selected = [ks_pts_not_selected; current_pt];
            else
                pts_not_selected = [pts_not_selected; current_pt];
            end
        else
            % Point was inside the 10% boundary, check if it was kept
            is_in_B = ismember(current_id, Set_B_IDs, 'rows');
            
            if is_in_B
                % 2. Selected by lineExclusion AND Kept by spikeMatch
                if is_ks
                    ks_pts_kept = [ks_pts_kept; current_pt];
                else
                    pts_kept = [pts_kept; current_pt];
                end
            else
                % 3. Selected by lineExclusion but REMOVED by spikeMatch
                if is_ks
                    ks_pts_redundant = [ks_pts_redundant; current_pt];
                else
                    pts_redundant = [pts_redundant; current_pt];
                end
            end
        end
    end
    
    % --- Start Plotting ---
    figure('Visible','on'); hold on; grid on; box on;
    h_legend = [];
    
    % Plot NOT selected points (Black 'x')
    if ~isempty(pts_not_selected)
        h_not_selected = scatter(pts_not_selected(:,1), pts_not_selected(:,2), 30, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'Not Selected (Outside Boundary)');
        h_legend(end+1) = h_not_selected;
    end
    
    % Plot REDUNDANT points (Green 'x')
    if ~isempty(pts_redundant)
        h_redundant = scatter(pts_redundant(:,1), pts_redundant(:,2), 40, [0 0.7 0], 'x', 'LineWidth', 2, 'DisplayName', 'Redundant (Removed by SpikeMatch)');
        h_legend(end+1) = h_redundant;
    end
    
    % Plot KEPT points (Green 'o')
    if ~isempty(pts_kept)
        h_selected = scatter(pts_kept(:,1), pts_kept(:,2), 40, [0 0.7 0], 'filled', 'DisplayName', 'Selected & Kept');
        h_legend(end+1) = h_selected;
    end
    
    % --- NEW KS PLOTTING (ON TOP) ---
    % Combine both 'x' types for KS
    ks_all_x = [ks_pts_not_selected; ks_pts_redundant];
    
    % Plot KS 'x' points (Blue 'x')
    if ~isempty(ks_all_x)
        h_ks_x = scatter(ks_all_x(:,1), ks_all_x(:,2), 40, 'b', 'x', 'LineWidth', 2, 'DisplayName', 'KS Coeff (Not Kept)');
        h_legend(end+1) = h_ks_x;
    end
    
    % Plot KS 'o' points (Blue 'o')
    if ~isempty(ks_pts_kept)
        h_ks_o = scatter(ks_pts_kept(:,1), ks_pts_kept(:,2), 40, 'b', 'filled', 'DisplayName', 'KS Coeff (Kept)');
        h_legend(end+1) = h_ks_o;
    end
    % --- END NEW KS PLOT ---
    
    % --- Plot 10% Trimmed Boundary Line ---
    keep_mask_10 = trim_struct_10.keep_mask;
    x_keep_10 = all_x(keep_mask_10);
    y_keep_10 = all_y(keep_mask_10);
    
    h_boundary = [];
    if ~isempty(x_keep_10)
        x_lim_Q2_10 = min(x_keep_10); y_lim_Q2_10 = max(y_keep_10);
        x_lim_Q4_10 = max(x_keep_10); y_lim_Q4_10 = min(y_keep_10);
    
        % Plot Q2 segment
        h_boundary = plot([x_lim_Q2_10, x_inter], [y_lim_Q2_10, y_inter], 'r-', 'LineWidth', 1.5, 'DisplayName', 'Boundary');
        % Plot Q4 segment
        plot([x_inter, x_lim_Q4_10], [y_inter, y_lim_Q4_10], 'r-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
    if ~isempty(h_boundary), h_legend(end+1) = h_boundary; end
    
    xL = xlim; yL = ylim;
    h_yline_inter = plot(xL, [y_inter y_inter], 'b:', 'LineWidth', 1.0, 'HandleVisibility','on', 'DisplayName', 'Intersection Line');
    h_xline_inter = plot([x_inter x_inter], yL, 'b:', 'LineWidth', 1.0, 'HandleVisibility','off'); 
    uistack([h_xline_inter, h_yline_inter], 'bottom'); 
    h_legend(end+1) = h_yline_inter;
    
    % Plot Intersection Point (should be on top)
    h_inter = plot(x_inter, y_inter, 'ko', 'MarkerFaceColor','y', 'MarkerSize',8, 'DisplayName','Intersection Point');
    h_legend(end+1) = h_inter;
    
    xlabel('Median Distance');
    ylabel('K-value');
    title(sprintf('Final Selection (Channel %s)', channelNumStr));
    legend(h_legend, 'Location','bestoutside');
    
    
    num_gauss_before = size(selected_gauss_10pct, 1);
    num_coeff_before = numel(unique(selected_gauss_10pct(:, 3)));
    
    num_gauss_after = size(select_spike_match, 1);
    num_coeff_after = numel(unique(select_spike_match(:, 3)));
    
    % --- NEW KS COUNT ---
    num_ks = numel(ks_coeff);
    
    % Format the text string
    txt_summary = sprintf([...
        '-- Before Redundancy (10%% Trim) --\n' ...
        '  Selected Gaussians: %d\n' ...
        '  Unique Coefficients: %d\n' ...
        '\n' ...
        '-- After Redundancy (SpikeMatch) --\n' ...
        '  Selected Gaussians: %d\n' ...
        '  Unique Coefficients: %d\n' ...
        '\n' ...
        '-- Special Sets --\n' ...                 % <-- NEW SECTION
        '  KS Coefficients: %d'], ...              % <-- NEW LINE
        num_gauss_before, num_coeff_before, ...
        num_gauss_after, num_coeff_after, ...
        num_ks); % Added num_ks
        
    % Position is [x, y, w, h] in normalized figure units.
    % We set a starting [x, y] and let 'FitBoxToText' handle the size.
    % --- MODIFIED X-POSITION from 0.83 to 0.85 ---
    annotation('textbox',[0.85, 0.15, 0.1, 0.1], 'String',txt_summary, ...
        'FitBoxToText','on', 'EdgeColor','k', 'BackgroundColor', 'w', ...
        'FontSize',9, 'VerticalAlignment','top', 'Margin', 5);
    
    
    filename3 = fullfile(folderName, sprintf('ch%s_10pct_final_selection.pdf', channelNumStr));
    exportgraphics(gcf, filename3, 'Resolution',300);
    close(gcf);
end