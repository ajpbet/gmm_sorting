function Idistcoeff_idx = plot_coeffs(pg,xg,wd_coeff,g,Idist,plotAll,medDist,ks_coeffs)

    % wd_coeff:     N x M matrix
    % g:            cell array of GMM models, length = M
    % Idist:        vector length M, normalized distance (will be unscaled for plotting)
    % plotAll:      boolean, true = plot all coeffs, false = only top 10 Idist

    % === PLOTTING MODE SELECTION ===
    % Change this variable to switch between plotting options:
    % 'idist'    - plot only Idist metric
    % 'meddist'  - plot only medDist metric  
    % 'both'     - plot both metrics
    plottingMode = 'both';  % Change this as needed
    
    num_coeff = size(wd_coeff,2);

    figure;
    histogram(Idist,15, 'Normalization','pdf');
    hold on;
    % x = linspace(min(Idist),max(Idist),200);
    % pd = fitdist(Idist,'Normal');
    [y,x] = ksdensity(Idist); 
   % y = pdf(pd, x);                  % evaluate the fitted PDF at each x
    
    plot(x, y, 'r-', 'LineWidth', 2);


    xlabel('Coefficient value');
    ylabel('Probability density');
    legend('Histogram', 'Normal fit PDF');
    title('Histogram with PDF overlay');
    grid on;    

    % Extract polyIdist from 64x2 cell array structure
    all_mu = cellfun(@(gm) gm.mu(:), g, 'UniformOutput', false);
    
    % Extract the polynomial indices from the second column
    polyIdist = medDist(:,2);  % This gives a cell array of single numbers
    polyIdist = cell2mat(polyIdist);  % Convert to a numeric vector
    
    mu_inputs = zeros(num_coeff,1);
    for k = 1:num_coeff
        if ~isempty(all_mu{k}) && polyIdist(k) <= length(all_mu{k})
            mu_inputs(k) = all_mu{k}(polyIdist(k));
        else
            mu_inputs(k) = NaN; % or some default
        end
    end
    % KDE + MSE (kde returns pdf and x-grid)
    [kde_pdf,kde_xf] = kde_est(wd_coeff,1:num_coeff);
    mse_vals = mean_square_error(pg,xg,kde_pdf,kde_xf,1:num_coeff);

    % choose which coefficients to plot
    if plotAll
        coeff_idx = 1:num_coeff;
        Idistcoeff_idx = 1:num_coeff;
    else
        [~,Idist_idx] = sort(Idist,"descend");
        Idistcoeff_idx = Idist_idx(1:min(14,num_coeff));  % top 10
        coeff_idx = unique([Idistcoeff_idx; ks_coeffs'], 'stable');
        %coeff_idx = Idistcoeff_idx;
    end

    % ----------- plotting loop -----------
    for ii = 1:length(coeff_idx)
        i = coeff_idx(ii);   % actual coefficient index
        
        if plotAll == false
            if ismember(i,Idistcoeff_idx)
                Idist_pt = "Idist mean ";
            else
                Idist_pt = "";
            end
              
            if ismember(i,ks_coeffs)
                ks_pt = "ks ";
            else
                ks_pt = "";
            end
        else
            Idist_pt = "";
            ks_pt = "";
        end
               % Create figure with proper layout from the start
        fig = figure('Name', sprintf('Coeff %d', i), 'NumberTitle', 'off', ...
                     'Position', [100, 100, 1200, 800]);  % Wider figure
        set(fig, 'WindowStyle', 'docked');
        
      % Create subplots with correct sizes from the beginning
        ax1 = subplot(2,3,[1 2]);  % Top row, left two-thirds
        ax2 = subplot(2,3,[4 5]);  % Bottom row, left two-thirds
        ax_table = subplot(2,3,[3 6]);  % Right column, spans both rows

        % --- subplot 1: histogram + KDE ---
        histogram(ax1, wd_coeff(:,i),1000,'Normalization','pdf');
        hold(ax1, 'on');
        plot(ax1, kde_xf{i},kde_pdf{i},'LineWidth',2);
        title(ax1, sprintf(Idist_pt+ks_pt+"Coeff %d — KDE (hist + kde)", i));
        xlabel(ax1, 'value'); ylabel(ax1, 'pdf');
        legend(ax1, 'hist','kde','Location','best');
        
        % >>> detect critical points on KDE (peaks & inflection) <<<
        xx_kde = xg{i}(:);
        yy_kde = pg{i}(:);
        yy_s = yy_kde; % keep as-is

        dy = gradient(yy_s, xx_kde);
        d2 = gradient(dy, xx_kde);

        % peaks: dy goes + -> - (positive to non-positive)
        pk_idx = find(dy(1:end-1) > 0 & dy(2:end) <= 0) + 1;
        % inflection: second derivative zero crossing (sign change)
        infl_idx = find(d2(1:end-1) > 0 & d2(2:end) <= 0) + 1;
        infl_idx2 = find(d2(1:end-1) < 0 & d2(2:end) >= 0) + 1;
        infl_idx = unique([infl_idx; infl_idx2]);

        % remove edge indices and very small peaks (threshold)
        eps_thresh = max(yy_s)*1e-3;
        pk_idx = pk_idx(pk_idx > 1 & pk_idx < length(xx_kde) & yy_s(pk_idx) > eps_thresh);
        infl_idx = infl_idx(infl_idx > 1 & infl_idx < length(xx_kde) & yy_s(infl_idx) > eps_thresh);

        % plot critical points on KDE
        if ~isempty(pk_idx)
            plot(ax1,xx_kde(pk_idx), yy_s(pk_idx), 'r^', 'MarkerFaceColor','r', 'MarkerSize',7, 'DisplayName','KDE peaks');
        end
        if ~isempty(infl_idx)
            plot(ax1,xx_kde(infl_idx), yy_s(infl_idx), 'gs', 'MarkerFaceColor','g', 'MarkerSize',7, 'DisplayName','KDE inflection');
        end
        hold(ax1,'off');

        % --- subplot 2: pdf + GMM components ---
        % show KDE overlay on PDF subplot (small offset so it's visible)
        plot(ax2,kde_xf{i},kde_pdf{i},'LineWidth',7,'Color',[0 0.4470 0.7410]);
        hold(ax2,'on')
        plot(ax2,xg{i},pg{i},'r','LineWidth',2);  
        hold(ax2,'on')

        gmm_model = g{i};
        mu = gmm_model.mu(:);
        sigma = sqrt(squeeze(gmm_model.Sigma(:)));
        w = gmm_model.ComponentProportion(:);

        % Normalize distance metrics for color mapping (per-coeff scaling)
        % Extract from first column of 64x2 cell array
        dm = medDist{i,1};
        if isempty(dm)
            dm = zeros(length(mu),1);
        end
        dm_norm = (dm - min(dm)) / (max(dm) - min(dm) + eps);
        cmap = parula(max(length(mu),64)); % ensure colormap has enough rows

        for k = 1:length(mu)
            comp_pdf = gmm_model.ComponentProportion(k) * normpdf(xg{i}, mu(k), sigma(k));
            c_idx = max(1, round(dm_norm(k) * (size(cmap,1)-1))+1);
            col = cmap(c_idx,:);
            plot(ax2,xg{i}, comp_pdf,'Color',col,'LineWidth',1.5);
            [~,pk_idx_comp] = max(comp_pdf);
            x_vals = xg{i};
            plot(ax2,x_vals(pk_idx_comp), comp_pdf(pk_idx_comp), 'o', 'MarkerFaceColor', col, ...
                 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
        end

        % colorbar for component distance metric (per-coeff)
        colormap(ax2, parula);
        caxis(ax2,[min(dm) max(dm)]);
        hcb = colorbar(ax2);
        hcb.Label.String = 'Distance Metric';

        % === PLOT METRICS BASED ON SELECTED MODE ===
        yl = ylim(ax2);
        range_mu = max(mu) - min(mu);
        
        % Get the polynomial indices for the highest values from second columns
        poly_idx_idist = medDist{i,2};    % polynomial index for highest Idist
        
        % Get the mu values for these specific polynomials
        mu_bar_idist = mu(poly_idx_idist);
        
        % Set different y-levels for arrows to avoid overlap
        y_arrow_idist = yl(2)*0.9;   % Higher level for Idist
        
        % Plot Idist metric
        if strcmpi(plottingMode, 'idist') || strcmpi(plottingMode, 'both')
            Idist_plot = Idist(i) * range_mu;   % unscaled for plotting
            dist_bar_idist = mu_bar_idist + Idist_plot;       % + distance
            dist_bar_minus_idist = mu_bar_idist - Idist_plot; % - distance

            % 3 vertical lines for Idist
            plot(ax2,[mu_bar_idist mu_bar_idist], yl, 'r--', 'LineWidth', 2, 'DisplayName', 'Idist μ');
            plot(ax2,[dist_bar_idist dist_bar_idist], yl, 'b--', 'LineWidth', 2, 'DisplayName', 'Idist bounds');
            plot(ax2,[dist_bar_minus_idist dist_bar_minus_idist], yl, 'b--', 'LineWidth', 2);
            
            % Arrow and label for Idist
            plot(ax2,[mu_bar_idist dist_bar_idist], [y_arrow_idist y_arrow_idist], '-^', 'Color', 'm', 'LineWidth', 2, ...
                'Marker', 'v', 'MarkerFaceColor', 'm', 'MarkerSize', 6, 'DisplayName', 'Idist distance');
            text(ax2,(mu_bar_idist+dist_bar_idist)/2, y_arrow_idist*1.02, sprintf('Idist: %.3g', Idist(i)),...
                'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color','m');
        end

        % >>> map KDE critical positions to the pdf grid (xg) <<<
        % find nearest xg indices for each KDE peak and inflection x
        isPeak_pdf = false(size(xg{i}));
        isInf_pdf  = false(size(xg{i}));
        % peaks
        for p = 1:numel(pk_idx)
            xp = xx_kde(pk_idx(p));
            [~, j] = min(abs(xg{i} - xp));
            isPeak_pdf(j) = true;
        end
        % inflections
        for p = 1:numel(infl_idx)
            xp = xx_kde(infl_idx(p));
            [~, j] = min(abs(xg{i} - xp));
            isInf_pdf(j) = true;
        end

        % compute metrics using pg (as you requested)
        pg_vec = pg{i};
        pgmax = max(pg_vec);
        Ipeak = 0; Iinf = 0;
        if pgmax > 0
            Ipeak = sum(pg_vec(isPeak_pdf)) / pgmax;
            Iinf  = sum(pg_vec(isInf_pdf)) / pgmax;
        end

        % plot markers on PDF too at matched x positions
        if any(isPeak_pdf)
            plot(ax2,xg{i}(isPeak_pdf), pg_vec(isPeak_pdf), 'r^','MarkerFaceColor','r','MarkerSize',7);
        end
        if any(isInf_pdf)
            plot(ax2,xg{i}(isInf_pdf), pg_vec(isInf_pdf), 'gs','MarkerFaceColor','g','MarkerSize',7);
        end

        % --- annotate MSE & Idist/medDist based on mode ---
        text(ax2,mean(xg{i})*0.65, max(pg{i})*1.05, sprintf('MSE = %.4f', mse_vals(i)),...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
        
        if strcmpi(plottingMode, 'idist') || strcmpi(plottingMode, 'both')
            text(ax2,mean(xg{i})*0.65, max(pg{i})*0.95, sprintf('Idist = %.4f', Idist(i)),...
                'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',[0 0.5 0]);
        end
        

        % >>> side summary box: counts and I metrics <<<
        cnt_peaks = sum(isPeak_pdf);
        cnt_infl  = sum(isInf_pdf);
        summary_str = {
            sprintf('Peaks: %d', cnt_peaks), ...
            sprintf('Inflections: %d', cnt_infl), ...
            sprintf('Ipeak = %.4f', Ipeak), ...
            sprintf('Iinf  = %.4f', Iinf) ...
            };
        % put textbox in upper-right of subplot 2
        xlims = xlim(ax2); ylims = ylim(ax2);
        tx = xlims(1) + 0.75*diff(xlims);
        ty = ylims(1) + 0.60*diff(ylims);
        text(ax2,tx, ty, summary_str, 'BackgroundColor','w','EdgeColor','k','FontSize',9);

        title(ax2, sprintf('Fig %d — PDF for coeff: %d (Mode: %s)', i, i, plottingMode));
        legend(ax2,'kde','pdf','','','','','','','location','best');
        hold(ax2,'off')
        
              % ----------- combined table: PolyIdx, Mean, Std, Weight, Idist, medDist -----------
        % ----------- Position table to the right of both subplots -----------
        idist_vals = medDist{i,1}; 
        
        info_str = cell(length(mu)+1,1);
        info_str{1} = sprintf('Poly |   Mean   |   Std   | Weight |  Idist   |');
        for k = 1:length(mu)
            info_str{k+1} = sprintf('%4d | %7.3f | %7.3f | %6.3f | %7.4f |', ...
                k, mu(k), sigma(k), w(k), idist_vals(k));
        end
        
        % Adjust font size based on number of components
        if length(mu) > 5
            font_size = 8;
        else
            font_size = 10;
        end
        
        text(ax_table, 0.05, 0.95, info_str, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontName', 'Courier', 'FontSize', font_size, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 2);
        axis(ax_table, 'off');

    end

    % ----------- summary tables -----------
    coeffs = coeff_idx(:);   % only the plotted ones
    table_summary_detect = table(coeffs, mse_vals(coeff_idx)', ...
        'VariableNames', {'coeffs','mseVals'});
    
    % Use the pre-computed highest values from Idist and Idist_med inputs
    table_summary_vals = table( ...
        coeffs,...
        mse_vals(coeff_idx)', ...                     
        mu_inputs(coeff_idx), ...                    
        Idist(coeff_idx), ...                        
        'VariableNames', {'coeffs','MSE','Mu','Idist'});
end

% --- helper functions ---
function [kde_pdf,kde_xf] = kde_est(wd_coeff,coeff_idx)
    kde_pdf = cell(1, length(coeff_idx));
    kde_xf  = cell(1, length(coeff_idx));
    for ii = 1:length(coeff_idx)
        i = coeff_idx(ii);
        [kde_pdf{ii},kde_xf{ii}] = kde(wd_coeff(:,i));
    end
end

function mse_vals = mean_square_error(pg,xg,kde_pdf,kde_xf,coeff_idx)
    mse_vals = zeros(1,length(coeff_idx));
    for ii = 1:length(coeff_idx)
        i = coeff_idx(ii);
        kde_pdf_interp = interp1(kde_xf{ii}, kde_pdf{ii}, xg{i}, 'linear');
        mse_vals(ii) = mean((normalize(kde_pdf_interp) - normalize(pg{i})).^2);
    end
    % Note: mse_vals returned is keyed to coeff_idx ordering
end