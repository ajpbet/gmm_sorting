function total_plots(pg,xg,wd_coeff,g,ks_out_full,ks_out,summary_table,spikes,all_ks,cluster_times,coeff_vector)
    %each norm gets one value
    %figures show how value interacts
    %gigure out how to remove gmm outlier
    %raw and tau spearman
    
    %first going to plot all kj and idist vs coeff
    % on each figure it will show a plot of k a plot of med dist
    % and a table of med idist ranked high to low with comp number
    % and similarly k val next to these ascending with associated comp number
    
    % graphs should show ranked plot of all k and idist
    % their comp values should also be present either on the x or form an 
    % annotation.

    %% plot only good coeffs or plot top coeffs
    plot_all = true;
    exc_combPdf = true;
    plot_none = false;

    %%
    % testing channels 332, 337,364,322
    lenKs = length(ks_out(:,1));
    coeff_nums = summary_table.("coeff num");
    ks_rankNums = flip(ks_out_full(:));
    ks_rankSel = ks_out(:,1);
    
    idistVals = summary_table.idist;
    max_k = summary_table.("max kj");
    
    [maxK_sort(:,2),maxK_sort(:,1)] = sort(max_k,"descend");
    
    [idist_sort(:,2),idist_sort(:,1)] = sort(idistVals,"descend");
    

    % change for idist over ks over lay ii double plot
    % need to lable
    figure;
    plot(maxK_sort(:,2))
    title(" max k vals per coeff")
    for k = 1:length(maxK_sort(:,2))
       text(k+0.1,maxK_sort(k,2)+0.1,num2str(maxK_sort(k,1))); 
    end

%% mathing top 3 vals
    maxLen = max(cellfun(@numel, summary_table.kj));
    kj_mat = cellfun(@(x) [x(:).' nan(1, maxLen - numel(x))], summary_table.kj, 'UniformOutput', false);
    kj_mat = cell2mat(kj_mat);

    [k_Mat,k_sortIdx]= sort(kj_mat, 2, 'descend');
    k_top3idx = k_sortIdx(:,1:3);

    medI_mat= cell2mat(cellfun(@(x) x(:).',summary_table.("med idist"),'UniformOutput',false));
    [medDist_sort,medDist_sortIdx] = sort(medI_mat,2,'descend');

%% exclude medDist values that are part of the combined pdf
    if exc_combPdf
        %combined_idx = cellfun(cellfun(@(x) x(:).',summary_table.M_comp,'UniformOutput',false));
        M_comb =cellfun(@(x) [x(:).' nan(1, maxLen - numel(x))], summary_table.M_comp, 'UniformOutput', false);
        combined_idx = cell2mat(M_comb);
        mask = true(size(medDist_sortIdx));  % initialize all ones
        numRows = numel(summary_table.M_comp);
        
        for r = 1:numRows
            idx = summary_table.M_comp{r};       % current row's indices to exclude
            mask(r, :) = ~ismember(medDist_sortIdx(r, :), idx);  % row-wise comparison
        end
        
        % Apply the mask
        medDist_sort = medDist_sort .* mask;
        medDist_sortIdx = medDist_sortIdx .* mask;
        
    end
%%
    count = 0;

    medDist_top3idx = zeros(size(medDist_sortIdx,1),3);
    for i = 1:size(medDist_sortIdx,1)
        for y = 1:size(medDist_sortIdx,2)
            if medDist_sortIdx(i,y) ~= 0 && ~isnan(medDist_sort(i,y))
               count = count+ 1;
               medDist_top3idx(i,count) = medDist_sortIdx(i,y);
               if count == 3
                   count = 0;
                   break
               end
            end
        end
    end
   % medDist_top3idx = medDist_sortIdx(:,1:min(3,length(medDist_sortIdx)));
    nRows = size(medDist_top3idx,1);
    comp_sharedIdx = zeros(nRows,1);  % result: one per row
    
    comp_shared = zeros(size(medDist_top3idx));

    for i = 1:nRows
        tempMed = medDist_top3idx(i,:);  % medDist top-3 for row i
        tempK   = k_top3idx(i,:);        % k top-3 for row i
        comp_shared(i,:) = tempMed.*ismember(tempMed,tempK);
        matchFound = false;
        for j = 1:length(tempMed)                       % loop over medDist in priority order
            if ismember(tempMed(j), tempK)
                comp_sharedIdx(i) = tempMed(j);  % store first match
                matchFound = true;
                break                        % move on to next row
            end
        end
        
        if ~matchFound
            comp_sharedIdx(i) = 0;           % no match at all
        end
    end

    %% use rankings for M comparison

    %% use weights of components

    %% rank and weights?


   %% first coeff comparison plot colors
    % first plot variables
    colors = cell(nRows,1);  % initialize 64x1, for example
    
    for i = 1:nRows
        matches_inTop3 = sum(ismember(medDist_top3idx(i,:), k_top3idx(i,:)));
        perfect_match = isequal(medDist_top3idx(i,:), k_top3idx(i,:));
    
        color = 'k';
        switch matches_inTop3
            case 0
                color = 'r';
            case 1
                color = 'y';
            case 2
                color = 'b';
            case 3
                if perfect_match
                    color = 'm';
                else
                    color = 'g';
                end
        end
        colors{i} = color;  % store it
    end
    
    % find which idist vals match a top 3 k val
     color_idistk = cell(nRows,1);
     idist_kmatch = zeros(nRows,1);
     ranking_Idist_match = zeros(nRows,1);
     for i = 1:nRows
        col = length(medDist_sort(i,:));
        for j = 1:col
            if ismember(medDist_sortIdx(i,j),k_top3idx(i,:))
                idist_kmatch(i) = medDist_sort(i,j);
                ranking_Idist_match(i) = j;
                break;
            end
        end
     end
     for i = 1:length(ranking_Idist_match)
         color_idistk{i} = 'r'; % default for idist loc >3
         switch ranking_Idist_match(i)
             case 1
                color_idistk{i} = 'm';
             case 2
                color_idistk{i} = 'g';
             case 3
                color_idistk{i} = 'b';
         end
     end
    

    %% coeff plots side by side
    figure;
    plot(idist_sort(:,2))
    title("idist vals per coeff")
    for k = 1:length(idist_sort(:,2))
       text(k+0.1,idist_sort(k,2)+0.1,num2str(idist_sort(k,1)),"Color",'b'); 
    end
    hold on;
    ks_d = flip(ks_out_full);
    ks_y = flip(all_ks);
    yyaxis right;
    plot(ks_y)
    hold on
    %% coeff scatter plot 1 - max values//color per matches
    
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
    xline(lenKs, '--', 'Color', 'k', 'LineWidth', 0.5, 'Alpha', 0.7);
    
    hold off;

        % scatter plot of both 
    figure;
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
    
    %% coefficient plot best idist/kj matching

    figure;
    hold on;
      
    % For each coefficient, plot IDist (x) vs KS (y)
    for k = 1:length(ks_y)
        coeff_num = ks_d(k);  % Coefficient number from KS ordering
        
        % Find the IDist value for this coefficient number
        idist_val = idist_kmatch(coeff_num);
        
        % Determine marker based on top 13 KS
        if k <= lenKs  % Since ks_y is sorted, first 13 are top
            % Top 13 KS - use dot
            marker = 'o';
            marker_face = color_idistk{coeff_num};
        else
            % Not top 13 KS - use x
            marker = 'x';
            marker_face = 'none';
        end
        
        % Plot the point
        scatter(idist_val, ks_y(k), 80, color_idistk{coeff_num}, marker, 'LineWidth', 1.5, ...
                'MarkerFaceColor', marker_face);
        
        % Label with coefficient number
        text(idist_val, ks_y(k)+0.01, num2str(coeff_num), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');
    end
    
    xlabel('IDist Values');
    ylabel('KS Values'); 
    title('IDist vs KS Values - highest idist with top 3 k)');
    grid on;
    
    % Match types (colors)
    h_match(1) = scatter(NaN, NaN, 80, 'm', 'o', 'filled', 'DisplayName', '1st mDist');
    h_match(2) = scatter(NaN, NaN, 80, 'g', 'o', 'filled', 'DisplayName', '2nd mDist');
    h_match(3) = scatter(NaN, NaN, 80, 'b', 'o', 'filled', 'DisplayName', '3rd mDist');
    h_match(4) = scatter(NaN, NaN, 80, 'r', 'o', 'filled', 'DisplayName', '>3 mDist');

    % KS selection (markers)
    h_ks(1) = scatter(NaN, NaN, 80, 'k', 'o', 'filled', 'DisplayName', 'KS selected');
    h_ks(2) = scatter(NaN, NaN, 80, 'k', 'x', 'LineWidth', 1.5, 'DisplayName', 'KS not selected');
    
    legend([h_match, h_ks], 'Location', 'best');
    hold off;



    %% code continues
    lenC = length(coeff_nums);
    %lenC = 20;
    %put in title ks ranking for each ranking rank by best ks values
    % is present in ks metric

    %find common overlap between k and idist
%% variable coefficient input
    if isempty(coeff_vector)
        coeff_vals = ks_rankSel;
    elseif isscalar(coeff_vector)
        switch coeff_vector
            case 1
                coeff_vals = ks_rankSel;
            case 2
                coeff_vals = idist_sort(1:length(ks_rankSel),1);
            case 3
                coeff_vals = idist_kmatch(1:length(ks_rankSel),1);
            case 4
                coeff_vals = ks_rankNums;
            case 5
                % test case
                coeff_vals = [1 2 3 4];
        end
    else
        coeff_vals = coeff_vector;
    end


    [kde_pdf,kde_xf] = kde_est(wd_coeff,1:lenC);
    mse_vals = mean_square_error(pg,xg,kde_pdf,kde_xf,1:lenC);
    
    %% main plots spikes etc.
    for i = 1:length(coeff_vals)
        % set up to go through coeffs by k rank
        coeff_num = coeff_vals(i);
        fig = figure('Name', sprintf('Coeff %d', coeff_num), 'NumberTitle', 'off', ...
                     'Position', [100, 100, 1200, 800]);
        set(fig, 'WindowStyle', 'docked');
        
        % Create main figure subplots
        ax1 = subplot(2,4,[1 2]);
        ax2 = subplot(2,4,[5 6]);
        ax_table = subplot(2,4,[3 4]);
        ax_table2 = subplot(2,4,[7 8]);
        
        % Your GMM model code
        gmm_model = g{coeff_num};
        mu = gmm_model.mu(:);
        sigma = sqrt(squeeze(gmm_model.Sigma(:)));
        w = gmm_model.ComponentProportion(:);
        dm = summary_table.("med idist"){coeff_num};
        
        % Initialize popup figure handle
        popupFig = [];
        usePopup = 0;
        % Check if we need a popup figure BEFORE the loop
        if nnz(comp_shared(coeff_num,:)) > 1 || plot_all
            popupFig = figure('Name', sprintf('Coeff %d extended spikes', coeff_num), ...
                             'NumberTitle', 'off', 'Position', [150, 150, 1000, 600]);
            % Create subplots for popup figure
            axn = cell(1, 3);
            axn{1} = subplot(1,3,1);
            axn{2} = subplot(1,3,2);
            axn{3} = subplot(1,3,3);
            usePopup = 1;
        end
        if plot_all == true
            comp_shared = medDist_top3idx;
        end
        if plot_none == true
            comp_shared = [0 0 0];
        end
        % Now start your main loop
        for p = 1:size(comp_shared,2)
            if comp_shared(coeff_num,p) ~= 0
                % Determine which axes to use for plotting
                if usePopup
                    % Use popup figure axes
                    current_ax = axn{p};
                    figure(popupFig); % Switch to popup figure
                else
                    % Use main figure axes
                    current_ax = ax1;
                    figure(fig); % Switch to main figure
                end
                
                axes(current_ax);
                hold on;
                
                mu_bestDist = mu(comp_shared(coeff_num,p));
                sigma_bestDist = sigma(comp_shared(coeff_num,p));
        
                range_st = mu_bestDist - 3*sigma_bestDist;
                range_end = mu_bestDist + 3*sigma_bestDist;
        
                spikeIdx = find(wd_coeff(:,coeff_num) >= range_st & wd_coeff(:,coeff_num) <= range_end);
                
                clus_colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];
                    [0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
                    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];
                    [0.586207 0.827586 0.310345];[0.965517 0.620690 0.862069];
                    [0.620690 0.758621 1.]];
                
                % Continue with the rest of your plotting code...
                isolated_spikes = spikes(spikeIdx,:);
                cluster_num = cluster_times(spikeIdx,1);
        
                [clusters_present,~,cls_idx] = unique(cluster_num);
                clust_counts = accumarray(cls_idx,1);
                percent_clust = (clust_counts/numel(cluster_num))*100;
            
                [sortedClsP,sortIdxCls] = sort(percent_clust,'descend');
                    
                top3_cluster = sortIdxCls(1:min(3,length(percent_clust)));
                top3_perc = sortedClsP(1:min(3,length(percent_clust)));
                
                otherPerc = sum(percent_clust) - sum(top3_perc);
        
                x = 1:size(isolated_spikes,2);
                mean_spikes = mean(isolated_spikes,1);
                std_spikes = std(isolated_spikes,0,1);
                select_spikes = length(isolated_spikes(:,1));
                total_spikes = length(spikes(:,1));
                perc_tot = 100*select_spikes/total_spikes;
                color_select = clus_colors(cluster_times(spikeIdx,1)+1,1:3);
        
                max_spikes_to_plot = 1000; % Adjust this number based on performance
                
                if size(isolated_spikes, 1) > max_spikes_to_plot
                    % Method 1: Random sampling (fastest)
                    selected_indices = randperm(size(isolated_spikes, 1), max_spikes_to_plot);
                    isolated_spikes_to_plot = isolated_spikes(selected_indices, :);
                    color_select_to_plot = color_select(selected_indices, :);
                    
                    % Update the counts for accurate percentages
                else
                    isolated_spikes_to_plot = isolated_spikes;
                    color_select_to_plot = color_select;
                end
                
                % Then in your plotting loop, use the subset:
                for k = 1:size(isolated_spikes_to_plot, 1)
                    plot(x, isolated_spikes_to_plot(k,:), 'Color', color_select_to_plot(k,:), 'LineWidth', 0.5);
                end
                
                pairStr=arrayfun(@(x,p) sprintf('%2d(%.2f%%)',x,p),top3_cluster,top3_perc,'UniformOutput',false);
                if otherPerc > 0
                    pairStr{end+1} = sprintf('Other(%.2f%%)', otherPerc);
                end
                perc_string =strjoin(pairStr,',');
        
                % Plot mean and std
                plot(x, mean_spikes, 'b', 'LineWidth', 2);
                plot(x, mean_spikes + std_spikes, 'r--', 'LineWidth', 1.5);
                plot(x, mean_spikes - std_spikes, 'r--', 'LineWidth', 1.5);
                title(sprintf('gauss %1d PofTotal: %.2f%% (%5d)',comp_shared(coeff_num,p), perc_tot,select_spikes)); % , %s',comp_shared(coeff_num,p), perc_tot,perc_string
                
                % Your legend code continues here...
                legend_labels = {};
                legend_handles = [];
                
                % Add top clusters to legend
                for i = 1:length(top3_cluster)
                    cluster_id = top3_cluster(i);
                    color_idx = cluster_id + 1;
                    if color_idx <= size(clus_colors, 1)
                        h = plot(NaN, NaN, 'Color', clus_colors(color_idx, :), 'LineWidth', 2);
                        legend_handles(end+1) = h;
                        legend_labels{end+1} = sprintf('clust %d (%.2f%%)', cluster_id, top3_perc(i));
                    end
                end
                
                % Add "Other" category
                if otherPerc > 0
                    other_clusters = setdiff(clusters_present, top3_cluster);
                    if ~isempty(other_clusters)
                        other_color = mean(clus_colors(other_clusters+1, :), 1);
                        h = plot(NaN, NaN, 'Color', other_color, 'LineWidth', 2);
                        legend_handles(end+1) = h;
                        legend_labels{end+1} = sprintf('Other (%.2f%%)', otherPerc);
                    end
                end
                
                % Add mean and std to legend
                h_mean = plot(NaN, NaN, 'b', 'LineWidth', 2);
                h_std = plot(NaN, NaN, 'r--', 'LineWidth', 2);
                legend_handles(end+1:end+2) = [h_mean, h_std];
                legend_labels{end+1} = 'Mean';
                legend_labels{end+1} = 'Mean ± STD';
                
                % Create the legend
                legend(legend_handles, legend_labels, 'Location', 'best');
                hold off;
            end
        end
        
        figure(fig);
        % >>> detect critical points on KDE (peaks & inflection) <<<
        xx_kde = xg{coeff_num}(:);
        yy_kde = pg{coeff_num}(:);
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

                % --- subplot 2: pdf + GMM components ---
        % show KDE overlay on PDF subplot (small offset so it's visible)
        plot(ax2,kde_xf{coeff_num},kde_pdf{coeff_num},'LineWidth',7,'Color',[0 0.4470 0.7410]);
        hold(ax2,'on')
        plot(ax2,xg{coeff_num},pg{coeff_num},'r','LineWidth',2);  
        hold(ax2,'on')

        if isempty(dm)
            dm = zeros(length(mu),1);
        end
        dm_norm = (dm - min(dm)) / (max(dm) - min(dm) + eps);
        cmap = parula(max(length(mu),64)); % ensure colormap has enough rows

        for k = 1:length(mu)
            comp_pdf = gmm_model.ComponentProportion(k) * normpdf(xg{coeff_num}, mu(k), sigma(k));
            c_idx = max(1, round(dm_norm(k) * (size(cmap,1)-1))+1);
            col = cmap(c_idx,:);
            plot(ax2,xg{coeff_num}, comp_pdf,'Color',col,'LineWidth',1.5);
            [~,pk_idx_comp] = max(comp_pdf);
            x_vals = xg{coeff_num};
            plot(ax2,x_vals(pk_idx_comp), comp_pdf(pk_idx_comp), 'o', 'MarkerFaceColor', col, ...
                 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
        end

        % colorbar for component distance metric (per-coeff)
        colormap(ax2, parula);
        caxis(ax2,[min(dm) max(dm)]);
        hcb = colorbar(ax2);
        hcb.Label.String = 'Distance Metric';

        isPeak_pdf = false(size(xg{coeff_num}));
        isInf_pdf  = false(size(xg{coeff_num}));
        % peaks
        for p = 1:numel(pk_idx)
            xp = xx_kde(pk_idx(p));
            [~, j] = min(abs(xg{coeff_num} - xp));
            isPeak_pdf(j) = true;
        end
        % inflections
        for p = 1:numel(infl_idx)
            xp = xx_kde(infl_idx(p));
            [~, j] = min(abs(xg{coeff_num} - xp));
            isInf_pdf(j) = true;
        end
        % compute metrics using pg (as you requested)
        pg_vec = pg{coeff_num};
        pgmax = max(pg_vec);
        Ipeak = 0; Iinf = 0;
        if pgmax > 0
            Ipeak = sum(pg_vec(isPeak_pdf)) / pgmax;
            Iinf  = sum(pg_vec(isInf_pdf)) / pgmax;
        end
        % plot markers on PDF too at matched x positions
        if any(isPeak_pdf)
            plot(ax2,xg{coeff_num}(isPeak_pdf), pg_vec(isPeak_pdf), 'r^','MarkerFaceColor','r','MarkerSize',7);
        end
        if any(isInf_pdf)
            plot(ax2,xg{coeff_num}(isInf_pdf), pg_vec(isInf_pdf), 'gs','MarkerFaceColor','g','MarkerSize',7);
        end

        text(ax2,mean(xg{coeff_num})*0.65, max(pg{coeff_num})*1.05, sprintf('MSE = %.4f', mse_vals(coeff_num)),...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
        
        title(ax2, sprintf('PDF for coeff: %d', coeff_num));
        legend(ax2,'kde','pdf','','','','','','','location','best');
        hold(ax2,'off')

        
        info_str = cell(length(mu)+1,1);
        info_str{1} = sprintf('Poly |   Mean   |   Std   | Weight |');
        for k = 1:length(mu)
            info_str{k+1} = sprintf('%4d | %7.3f | %7.3f | %6.3f |', ...
                k, mu(k), sigma(k), w(k));
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

        info_str2 = cell(length(k_Mat(1,:))+1,1);
        info_str2{1} = sprintf('K_Values|K_Values_IDX| MedDist | MedI_IDX |');
        for k = 1:length(k_Mat(1,:))
            info_str2{k+1} = sprintf('%7.3f | %7d | %7.3f | %7d |', ...
            k_Mat(coeff_num,k), k_sortIdx(coeff_num,k), medDist_sort(coeff_num,k), medDist_sortIdx(coeff_num,k));
        end

        text(ax_table2, 0.05, 0.95, info_str2, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontName', 'Courier', 'FontSize', font_size, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 2);
        axis(ax_table2, 'off');
        
    end
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