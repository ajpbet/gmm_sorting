function [select_all, ks_coeff,select_spike_match] = GMM_extract(spikes, par, filename_mat,basename)
    % If GMM .mat file doesn't exist, generate and save it
    if ~exist(filename_mat, 'file')
        [inspk, ks_coeff, ks_out, ks_out_full, all_ks] = wave_features(spikes, par);
        [g, xg, pg, coeffs, medDist, summary_table] = GMM_full(spikes, par);

        save(filename_mat, "summary_table", "g", "pg", "xg", ...
             "coeffs", "inspk", "ks_coeff", "ks_out", ...
             "ks_out_full", "all_ks");
    else
        load(filename_mat);
        [inspk, ks_coeff, ks_out, ks_out_full, all_ks] = wave_features(spikes, par);
        save(filename_mat, "summary_table", "g", "pg", "xg", ...
             "coeffs", "inspk", "ks_coeff", "ks_out", ...
             "ks_out_full", "all_ks");
    end

    %
    folderName = 'results';
    folderPlots = fullfile(folderName, 'plots');
    folderSpikeMatch = fullfile(folderName,'spikeMatch');

    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    if ~exist(folderPlots, 'dir')
        mkdir(folderPlots);
    end
    if ~exist(folderSpikeMatch, 'dir')
        mkdir(folderSpikeMatch);
    end    

        %% filter GMM results by weight (K, M_comp, medDist, polyID_M)
    g_init = cell(1,length(g));
    medDist_cell = cell(1,length(g));
    kj_mat = cell(1,length(g));
    M_comp = cell(1,length(g));
    polyID = cell(1,length(g));
    medD_sortAll = cell(1,length(g));
    medD_sIdxAll = cell(1,length(g));
    kj_sortAll = cell(1,length(g));
    kj_sIdxAll = cell(1,length(g));

    idistVals = zeros(length(g),1);
    max_k = zeros(length(g),1);

    medDist_init = summary_table.("med idist");
    Mcomp_init = summary_table.M_comp;
    kj_init = summary_table.kj;
    g_init = g;

    for i = 1:length(g)
        g_in = g{i};
        mu = g_in.mu(:);                       
        s  = sqrt(squeeze(g_in.Sigma(:)));      
        a  = g_in.ComponentProportion(:);
        polyID_in = (1:numel(mu))';   % assign original indices from g
      
        keep = a > 0.005;   %way to make sure it is not isolated % 2.5% cutoff
        mu = mu(keep);
        s  = s(keep);
        a  = a(keep);
        polyID{i} = polyID_in(keep);
        Knew = numel(mu);               % actual surviving components
        a_new = a / sum(a);
        s_new = reshape(s.^2, 1, 1, Knew);
        
        medDist_cell{i} = medDist_init{i}(keep);
        kj_mat{i} = kj_init{i}(keep);
        g{i} = gmdistribution(mu,s_new,a_new);
        M_comp{i} = Mcomp_init{i}(ismember(Mcomp_init{i},polyID{i}));
        pg2 = pdf(g{i}, xg{i});

        [medD_sortAll{i},medD_sIdxAll{i}] = sort(medDist_init{i},'Descend');
        [kj_sortAll{i},kj_sIdxAll{i}] = sort(kj_init{i},'Descend');
        
        if Knew > 0
            idistVals(i) = max(medDist_cell{i});
            max_k(i) = max(kj_mat{i});
        end

    end

        exc_combPdf = true;
    plot_none = false;

    polyID_M = polyID;
    %%
    % testing channels 332, 337,364,322
    lenKs = length(ks_out(:,1));
    coeff_nums = summary_table.("coeff num");
    ks_rankNums = flip(ks_out_full(:));
    ks_rankSel = ks_out(:,1);
    
    % fix 

    
    [maxK_sort(:,2),maxK_sort(:,1)] = sort(max_k,"descend");
    
    [idist_sort(:,2),idist_sort(:,1)] = sort(idistVals,"descend");
    

    % change for idist over ks over lay ii double plot
    % need to lable
    % fig_k = figure;
    % plot(maxK_sort(:,2))
    % title(" max k vals per coeff")
    % for k = 1:length(maxK_sort(:,2))
    %    text(k+0.1,maxK_sort(k,2)+0.1,num2str(maxK_sort(k,1))); 
    % end
    % filename_kv = fullfile(folderName,sprintf('ch%s_maxKvalsPerCoef.png', channelNum));
  %  exportgraphics(fig_k, filename_kv, 'Resolution', 300);

%% macthing top 3 vals
    maxLen = max(cellfun(@numel, kj_mat));
    
    % Sort each numeric array in descending order
    k_sort = cell(size(kj_mat));
    % Map sorted indices to polyIDs
    k_sortIdx = cell(size(kj_mat));
    for r = 1:numel(k_sort)
        vals = kj_mat{r};
        ids  = polyID_M{r};
        
        % Exclude NaNs
        mask = ~isnan(vals);
        vals_valid = vals(mask);
        ids_valid  = ids(mask);
        
        [vals_sorted, idx] = sort(vals_valid, 'descend');  % sort original numeric values
        k_sort{r} = vals_sorted;
        k_sortIdx{r} = ids_valid(idx);        % map sorted indices to polyIDs
    end
    
    % Top 3 polyIDs for plotting / comparison
    k_top3idx = zeros(length(k_sortIdx),3);
    for r = 1:length(k_sortIdx)
        nTop = min(3,length(k_sortIdx{r}));
        k_top3idx(r,1:nTop) = k_sortIdx{r}(1:nTop);
    end

    medDist_sort = cell(size(medDist_cell));
    medDist_sortIdx = cell(size(medDist_cell));
    
    for r = 1:numel(medDist_cell)
        vals = medDist_cell{r};
        ids  = polyID_M{r};
        
        % Exclude NaNs
        mask = ~isnan(vals);
        vals_valid = vals(mask);
        ids_valid  = ids(mask);
        
        % Sort the valid values
        [vals_sorted, idx] = sort(vals_valid, 'descend');
        
        medDist_sort{r} = vals_sorted;          % store sorted medDist values
        medDist_sortIdx{r} = ids_valid(idx);    % map sorted indices to polyIDs
    end
%% exclude medDist values that are part of the combined pdf
    numRows = numel(M_comp);
    
    medDist_sort_masked    = cell(size(medDist_sort));
    medDist_sortIdx_masked = cell(size(medDist_sortIdx));
    
    k_sort_masked = cell(size(k_sort));
    k_sortIdx_masked = cell(size(k_sortIdx));

    for r = 1:numRows
        vals = medDist_sort{r};
        ids  = medDist_sortIdx{r};
        idx_exclude = M_comp{r};  % indices to exclude
        
        mask = ~ismember(ids, idx_exclude) & ~isnan(vals);
        
        medDist_sort_masked{r}    = vals(mask);
        medDist_sortIdx_masked{r} = ids(mask);
    end
    
    for r = 1:numRows
        vals = k_sort{r};
        ids  = k_sortIdx{r};
        idx_exclude = M_comp{r};  % indices to exclude
        
        mask = ~ismember(ids, idx_exclude) & ~isnan(vals);
        
        k_sort_masked{r}    = vals(mask);
        k_sortIdx_masked{r} = ids(mask);
    end
    % Replace originals
    k_sort = k_sort_masked;
    k_sortIdx = k_sortIdx_masked;
    medDist_sort    = medDist_sort_masked;
    medDist_sortIdx = medDist_sortIdx_masked;


    %% kv knee(all present gauss after combined gaussian removed)
    % medDist_vec has struc [medDist,coeff#,component(gauss)#]
    [k_select,k_lSel,kDist_vec] = processKvKnee(k_sort, k_sortIdx);

    % process k knee when gaussians not near a peak are exluded
 %   [k_sel_NoPk,k_lSel_NoPk,kDist_vec_NoPk] = KvKneeNoPkExc(k_sort, k_sortIdx,pg,g_init,xg);
    %pltSelect(k_lSel,kDist_vec,'61',folderPlots,1);


        %% meddist knee (all present gauss after combined gaussian removed)
    [medDist_select, medDist_lSel,medDist_vec] = processMedDistKnee(medDist_sort, medDist_sortIdx);    %% kv plot all vals with knee (all present gauss after combined gaussian removed)
    % medDist_vec has struc [medDist,coeff#,component(gauss)#]
    % original polyIds preserved in gauss(id)
   % pltSelect(medDist_lSel,medDist_vec,'61',folderPlots,3);

    % process medDist when gaussians not near a peak are excluded
    % [medD_sel_noPk,medD_lSel_noPk,medD_vec_noPk] = medDistKneeNoPkExc(medDist_sort, medDist_sortIdx, ...
    % pg,g_init,xg);
%    pltDualKneePlot(k_lSel, kDist_vec, medDist_lSel, medDist_vec, 'KDist', 'medDist', basename, folderPlots)

       %% Exclusion criterion
   % with all peaks
    plot_and_save = true; % Or false
    
    [select_gauss,select_gauss1pct,select_gauss2p5pct,...
        select_gauss_5pct, select_gauss_10pct, select_all, ...
     x_vals, y_vals, gauss_ids, x_inter, y_inter] = lineExclusion(medDist_vec, ...
        medDist_select, k_select, kDist_vec, folderPlots, basename);
    
    % Get the boundary lines for the 10% trim
    threshQ2_10pct = select_all.trim10_0pct.threshQ2;
    threshQ4_10pct = select_all.trim10_0pct.threshQ4;
    
    select_spike_match = spikeMatch(select_gauss_10pct, ...
        threshQ2_10pct, threshQ4_10pct, g_init, coeffs, ...
        folderSpikeMatch, basename, 5, 0.9);
    pltDualKneePlot(k_lSel, kDist_vec, medDist_lSel, medDist_vec,select_gauss_10pct,select_spike_match, 'KDist', 'medDist', basename, folderPlots)


    % This plot now shows the final combined results
    plotSelectedVsAll(x_vals, y_vals, gauss_ids, ...
        x_inter, y_inter, ...
        select_all.trim10_0pct, select_gauss_10pct, ...
        select_spike_match, ...
        folderSpikeMatch, basename, plot_and_save,ks_coeff);

    
    %     k_select, folderPlots, basename);

   % without variables not near a peak

    % [select_gaussNoPk,select_gauss1pctNoPk,select_gauss2p5pctNoPk,...
    %     select_gauss_5pctNoPk, select_gauss_10pctNoPk,select_allNoPk] = ...
    %     lineExclusion(medD_vec_noPk,medD_sel_noPk, ...
    %     k_sel_NoPk,kDist_vec_NoPk);

    
end