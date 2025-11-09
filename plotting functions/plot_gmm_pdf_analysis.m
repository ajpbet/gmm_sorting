function [summary_list] = plot_gmm_pdf_analysis(gmm_cells, pdf_y_cells, dist_indices, x_range_cells, highlight_comp_idx)
% PLOT_GMM_PDF_ANALYSIS plots multiple indexed GMMs and their components 
% from cell arrays, assuming the overall PDF is provided as a data array.
%
% This function is designed for 1-dimensional GMM analysis and is compatible
% with MATLAB's gmdistribution object property names.
%
% Inputs:
%   gmm_cells           - Cell array (e.g., 1x64) where each cell contains a gmdistribution 
%                         object.
%   pdf_y_cells         - Cell array (e.g., 1x64) where each cell contains the 
%                         overall PDF data as a numerical Y-vector. (MANDATORY Y-DATA)
%   dist_indices        - An ARRAY of scalar indices (e.g., [1, 5, 10]) specifying 
%                         which distributions to select and plot.
%   x_range_cells       - Cell array (e.g., 1x64) where each cell contains the 
%                         x-axis domain vector (xg). This vector MUST correspond 
%                         to the data in pdf_y_cells. (MANDATORY X-DATA)
%   highlight_comp_idx  - OPTIONAL: An array of component indices (e.g., [11, 12]) 
%                         to be highlighted in red. The component index must be 
%                         <= K (number of components).
%
% Output:
%   summary_list  - A cell array where each cell contains the summary 
%                   structure for the corresponding plot index.
%
% A figure is automatically generated and saved for each index in dist_indices.
%
    if nargin < 4 || isempty(x_range_cells)
        error('x_range_cells (the X-data for your PDF) is now a mandatory input.');
    end
    
    % Handle optional input for the component to highlight
    if nargin < 5
        highlight_comp_idx = []; % Default to no component highlighted
    end

    num_plots = length(dist_indices);
    summary_list = cell(1, num_plots); % Initialize output cell array for all summaries
    disp(['Starting analysis for ', num2str(num_plots), ' distributions.']);
    
    % Define the new red color for the highlighted component
    highlight_color = [0.8 0 0.2]; % Same as the peak color for consistency

    % Loop through every requested distribution index
    for k = 1:num_plots
        
        dist_idx = dist_indices(k); % Get the current original index
        
        disp(['--- Processing Distribution Index: ', num2str(dist_idx), ' ---']);
        % Check if index is valid
        if dist_idx > length(gmm_cells) || dist_idx <= 0
            warning(['Index ', num2str(dist_idx), ' is out of bounds. Skipping.']);
            continue;
        end
        
        % --- 1. Input Selection and Data Extraction ---
        
        gmm = gmm_cells{dist_idx};
        
        % CRITICAL: Direct assignment of X and Y data vectors
        y_pdf = pdf_y_cells{dist_idx};
        x_range = x_range_cells{dist_idx};
        
        % Check data vector dimensions
        if length(y_pdf) ~= length(x_range)
             error(['PDF X and Y data vectors for index ', num2str(dist_idx), ' must have the same length.']);
        end
        
        % Extract GMM properties (gmdistribution compatible)
        weights = gmm.ComponentProportion(:);
        K = length(weights);
        means = gmm.mu(:);
        
        % Handle Variance Extraction (Sigma 1x1xK format)
        if ndims(gmm.Sigma) == 3 && size(gmm.Sigma, 1) == 1 && size(gmm.Sigma, 2) == 1
            variances = squeeze(gmm.Sigma);
        elseif isvector(gmm.Sigma) && length(gmm.Sigma) == K
            variances = gmm.Sigma(:);
        else
            error(['GMM Sigma format not supported for 1D analysis at index ', num2str(dist_idx), '.']);
        end
        stds = sqrt(variances); % Standard deviations
        
        
        figure('Color', 'w'); % Create a new figure for this plot
        hold on;
        
        % --- 2. Calculate Summary Statistics (The "Summary Table" Logic) ---
        
        % Find the peak of the input PDF (using the Y-data array)
        [y_peak_val, peak_idx] = max(y_pdf); % Capture both value and index
        x_peak = x_range(peak_idx);
        
        % Determine the 'width' W: 2 * standard deviation of the DOMINANT GMM component
        [~, dominant_idx] = max(weights);
        W = 2 * stds(dominant_idx); 
        
        x_start = x_peak - W / 2;
        x_end = x_peak + W / 2;
        
        % M_comp: Indices of Gaussians whose mean lies within the calculated width
        M_comp = find(means >= x_start & means <= x_end);
        
        % Create the temporary output summary structure
        temp_summary.Index = dist_idx;
        temp_summary.M_comp = M_comp;
        temp_summary.x_peak = x_peak;
        temp_summary.y_peak = y_peak_val;
        temp_summary.Width = W;
        temp_summary.x_start = x_start;
        temp_summary.x_end = x_end;
        
        % --- 3. Plotting Setup ---
        
        % Define colors
        main_pdf_color = [0 0 0];
        middle_comp_color = [0 0.5 0.8]; % Blue
        side_comp_color = [0.8 0.5 0];  % Orange
        shade_color = [0.9 0.9 0.9];    % Light Gray0
        
        max_y = max(y_pdf) * 1.1; 
        
        % --- 4. Plot Overall PDF and Shaded Area ---
        
        % Shade the "middle gaussian" width area
        fill_x = [x_start, x_end, x_end, x_start];
        fill_y = [0, 0, max_y, max_y];
        patch(fill_x, fill_y, shade_color, 'EdgeColor', 'none', ...
              'FaceAlpha', 0.5, 'DisplayName', 'Middle Region Width');
        % Plot the overall PDF (using the actual X and Y data)
        plot(x_range, y_pdf, 'LineWidth', 2.5, 'Color', main_pdf_color, ...
             'DisplayName', 'Overall PDF');
        
        % --- 5. Plot Individual Gaussian Components and Analysis Lines ---
        
        is_highlight_plotted = false; % Flag for the new legend entry

        for i = 1:K
            mu_i = means(i);
            sigma_i = stds(i);
            weight_i = weights(i);
            
            % Calculate component PDF against the user-provided X-range
            y_comp = weight_i * normpdf(x_range, mu_i, sigma_i);
            
            % Determine line style (Dotted if weight < 0.05)
            if weight_i < 0.005
                line_style = ':';
                line_width = 1;
                weight_label = ' (Weight < 0.05 - Dotted)';
            else
                line_style = '-';
                line_width = 1.5;
                weight_label = ' (Weight \geq 0.05 - Solid)';
            end
            
            % Determine color: Red if highlighted, Blue if in M_comp, Orange otherwise
            if i == highlight_comp_idx(k)
                line_color = highlight_color;
                group_label = ' (HIGHLIGHTED)';
                is_highlight_plotted = true; % Set flag
                plot([mu_i, mu_i], [0, max_y], ...
                     'Color', highlight_color, 'LineStyle', '--', 'LineWidth', 1.5, ...
                     'DisplayName', ['Comp ' num2str(i) ' Mean (x=' num2str(mu_i, '%.3f') ')']);
            elseif ismember(i, M_comp)
                line_color = middle_comp_color;
                group_label = ' (Middle Group)';
            else
                line_color = side_comp_color;
                group_label = ' (Side Group)';
            end
            
            % Plot the component
            plot(x_range, y_comp, line_style, 'Color', line_color, ...
                 'LineWidth', line_width, ...
                 'DisplayName', ['Comp ' num2str(i) group_label weight_label]);
            
            % Plot mean line for clarity
            % plot([mu_i, mu_i], [0, max(y_comp)], ...
            %      'Color', line_color, 'LineStyle', line_style, 'LineWidth', 0.5);
        end
        
        % --- 6. Plot Peak and Width Boundaries ---
        
        % Plot the PDF Peak as a dot (the required change)
        plot(x_peak, y_peak_val, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0 0.2], ...
             'MarkerEdgeColor', [0.8 0 0.2], 'LineStyle', 'none', 'DisplayName', 'PDF Peak');
        % Plot Width Boundaries
        plot([x_start, x_start], [0, max_y], 'Color', [0 0 0], 'LineWidth', 1, ...
             'LineStyle', '-', 'DisplayName', 'Width Boundary');
        plot([x_end, x_end], [0, max_y], 'Color', [0 0 0], 'LineWidth', 1, ...
             'LineStyle', '-', 'HandleVisibility', 'off'); % Hide from legend
         
        % --- 7. Final Polish and Saving ---
        
        title(['GMM Component Analysis for Distribution Index ', num2str(dist_idx)]);
        xlabel('X');
        ylabel('Probability Density');
        
        % --- Manual Legend using invisible placeholders ---
        hold on;
        h_pdf     = plot(nan, nan, 'k-',  'LineWidth', 2.5);            % PDF (black solid)
        h_width   = plot(nan, nan, 'k-',  'LineWidth', 1);              % Width (black thin)
        h_peak    = plot(nan, nan, 'o',   'MarkerSize', 8, ...
                         'MarkerFaceColor', [0.8 0 0.2], ...
                         'MarkerEdgeColor', [0.8 0 0.2]);               % Peak (red circle)
        h_highlight = plot(nan, nan, '-', 'Color', highlight_color, ... % NEW: Highlighted Component
                           'LineWidth', 1.5);
        
        h_main    = plot(nan, nan, '-',   'Color', [0 0.5 0.8], ...
                         'LineWidth', 1.5);                             % Components (Main Mixture)
        h_not_main    = plot(nan, nan, '-',   'Color', [0.8 0.5 0], ...
                         'LineWidth', 1.5);
        h_removed = plot(nan, nan, ':',   'Color', [0.8 0.5 0], ...
                         'LineWidth', 1.5);                             % Removed by Weight (dotted orange)
        
        % Conditionally include the highlighted component in the legend
        if is_highlight_plotted
            legend_handles = [h_pdf, h_width, h_peak, h_highlight,h_not_main, h_main, h_removed];
            legend_strings = {'PDF', 'Width', 'Peak', 'Highlighted Gaussian', 'Gaussians','Gaussians (Main Mixture)', 'Removed by Weight'};
        else
            % If no component was highlighted, use the original legend structure
            legend_handles = [h_pdf, h_width, h_peak, h_main, h_removed];
            legend_strings = {'PDF', 'Width', 'Peak', 'Components (Main Mixture)', 'Removed by Weight'};
        end
        
        legend(legend_handles, legend_strings, ...
               'Location', 'NorthEastOutside', 'Interpreter', 'none');

        grid on;
        ylim([0, max_y]);
        xlim([x_range(1), x_range(end)]);
        hold off;
        
        % Save the figure using the index in the filename
        filename = ['GMM_PDF_Analysis_Plot_', num2str(dist_idx), '.png'];
        saveas(gcf, filename);
      %  close(gcf); % Close the figure after saving
        
        % Store the summary structure in the list
        summary_list{k} = temp_summary;
        
        disp(['Figure successfully saved as: ' filename]);
    end
    
    disp('--- Batch analysis complete. Summary list returned. ---');
end