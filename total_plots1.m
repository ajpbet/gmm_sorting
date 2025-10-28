function total_plots(pg,xg,wd_coeff,g,ks_out_full,ks_out,summary_table,spikes,all_ks,cluster_times,coeff_vector, inputFile)
  % graphs should show ranked plot of all k and idist
    % their comp values should also be present either on the x or form an 
    % annotation.
    chStr = regexp(inputFile, '\d+', 'match');  % finds all numbers
    channelNum = chStr{2};  
    folderName = sprintf('ch%s', channelNum);
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end


    % filter GMM results by weight (K, M_comp, medDist, polyID_M)
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

    
        %% kv distributions with and without combined
    [maxK_sort(:,2),maxK_sort(:,1)] = sort(max_k,"descend");
    [kv_out,polyID_kv] = kv_dist(summary_table,channelNum);

        % filtering out which kvals to use no low weight or mcomp

                %% plot k values in a table with all coefficeints
    thresh_kv = 0;
    tab_gauss(summary_table,channelNum);


    %% plot max kv per coeff
    kv_per_coeff(ks_out,folderName,channelNum,maxK_sort);


    %% match top 3 vals of k with highest medDist
    [k_top3idx, medDist_top3idx, comp_shared, comp_sharedIdx] = processTop3Values(kj_mat, polyID_M, ...
      medDist_cell, M_comp, true);

    %% color control for plotting
    [idist_kmatch,color_idistk,poly_match_idistK,ranking_Idist_match] = early_color_plot_proc(medDist_top3idx,k_top3idx);

    %% idist selection knee
    [coeff_idist,idist_select] = idist_knee_select(idist_kmatch);

    idist_kmatch_lSel = length(idist_select(:,1));

    %% side by side coeff plot for idist and ks vals
    [ks_d,ks_y,r_sorted_idist,r_ind_idist] = plot_idist_ks(sorted_idist, ind_idist, idist_kmatch_lSel, ks_out_full, ...
        all_ks, lenKs, folderName, channelNum);

    %% scatter plot of max idist v. ks
     scatter_max_idist_ks(ks_d,ks_y,lenKs)

end