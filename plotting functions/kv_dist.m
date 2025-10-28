function kv_dist(summary_table, channelNum, folderName,g)
    % Extract kj values from the summary table (64x1 cell array, each with 8 values)
    kj_values = summary_table.kj;
    
    % Initialize arrays for kj values within and outside M_comp
    kj_within = [];
    kj_outside = [];
    
    % For each row, check which kj values are within that row's M_comp
    for i = 1:height(summary_table)
        current_kj = kj_values{i};  % Get the kj values from this row
        current_M_comp = summary_table.M_comp{i};  % Get the M_comp indices for this row

        g_in = g{i};
        comp_proportions = g_in.ComponentProportion;
        keep = comp_proportions > 0.005;
        
        
        % Apply the filter
        filtered_kj = current_kj(keep);
        filtered_positions = 1:length(current_kj);
        filtered_positions = filtered_positions(keep);
        
        % For this row's filtered kj values, check which positions are in this row's M_comp
        is_within_current = ismember(filtered_positions, current_M_comp);
        
        % Separate the filtered kj values for this row
        kj_within = [kj_within; filtered_kj(is_within_current)];
        kj_outside = [kj_outside; filtered_kj(~is_within_current)];
    end
    
    % Create common range for KDE comparison
    all_kj = [kj_within; kj_outside];
    if isempty(all_kj)
        warning('No kj values found for channel %s', channelNum);
        return;
    end
    
    min_val = min(all_kj);
    max_val = max(all_kj);
    common_range = linspace(min_val, max_val, 1000);
    
    % Calculate KDEs for both groups
    if ~isempty(kj_within)
        kde_within = ksdensity(kj_within, common_range);
    else
        kde_within = zeros(size(common_range));
    end
    
    if ~isempty(kj_outside)
        kde_outside = ksdensity(kj_outside, common_range);
    else
        kde_outside = zeros(size(common_range));
    end
    
    % Create figure
    fig = figure('Position', [100, 100, 800, 600]);
    hold on;
    
    % Plot KDEs
    if ~isempty(kj_within)
        plot(common_range, kde_within, 'b-', 'LineWidth', 2, 'DisplayName', 'Within M\_comp');
    end
    
    if ~isempty(kj_outside)
        plot(common_range, kde_outside, 'r-', 'LineWidth', 2, 'DisplayName', 'Outside M\_comp');
    end
    
    % Customize plot
    xlabel('kj Values', 'FontSize', 12);
    ylabel('Probability Density', 'FontSize', 12);
    title(sprintf('Distribution of kj Values - Channel %s', channelNum), 'FontSize', 14);
    legend('show', 'Location', 'best');
    grid on;
    
    % % Add some statistics to the plot
    % stats_text = sprintf('Within M\\_comp: %d values\nOutside M\\_comp: %d values\n(Components with proportion < 0.005 excluded)', ...
    %     length(kj_within), length(kj_outside));
    % 
    % text(0.05, 0.95, stats_text, ...
    %     'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white', ...
    %     'VerticalAlignment', 'top');
    % 
    % Create directory if it doesn't exist
    % Create absolute path
    absolute_folder = fullfile(pwd, folderName);
    
    % Create directory if it doesn't exist
    if ~exist(absolute_folder, 'dir')
        mkdir(absolute_folder);
    end
    
    % Save with absolute path
    saveas(fig, fullfile(absolute_folder, sprintf('kj_distributions_ch%s.png', channelNum)));
    hold off;
    
    % Display summary statistics
    fprintf('Channel %s Summary:\n', channelNum);
    fprintf('  Total kj values: %d\n', length(all_kj));
    if ~isempty(kj_within)
        fprintf('  Within M_comp:  %d values, mean=%.4f, std=%.4f\n', ...
            length(kj_within), mean(kj_within), std(kj_within));
    else
        fprintf('  Within M_comp:  0 values\n');
    end
    if ~isempty(kj_outside)
        fprintf('  Outside M_comp: %d values, mean=%.4f, std=%.4f\n', ...
            length(kj_outside), mean(kj_outside), std(kj_outside));
    else
        fprintf('  Outside M_comp: 0 values\n');
    end
    fprintf('  (Components with proportion < 0.005 excluded)\n');
end