function [g,g_filt,xg,pg,coeffs,medDist,best_values,summary_out] = GMM_full(spikes,par)
    coeffs = haar_coeffs(spikes,par);
 

    m = size(coeffs,2);

   
    medDist = cell(2,m);

    g      = cell(1,m); 
    g_filt = cell(1,m);         
    xg     = cell(1,m);           
    pg     = cell(1,m);           

    summary_gmm = cell(21,m);

    coeff_interest = 1;
    for j = 1:m
        Mj = GMM_basic1D(coeffs(:,j),par);

        [kj_max,kj_idx] = sort(Mj.dist_kj,'descend');
        [idist_max,idist_idx] = sort(Mj.Dij,'descend');


        summary_gmm{1,j} = j;
        summary_gmm{2,j} = Mj.dist_kj;
        summary_gmm{3,j} = kj_max(1);
        summary_gmm{4,j} = median(Mj.dist_kj);
        summary_gmm{5,j} = kj_idx(1:2);
        summary_gmm{6,j} = Mj.Dij;
        summary_gmm{7,j} = Mj.Idist;
        summary_gmm{8,j} = Mj.med_dist;
        summary_gmm{9,j} = Mj.mu_mpdf;
        summary_gmm{10,j} = Mj.s_mpdf;
        summary_gmm{11,j} = Mj.Kpeaks;
        summary_gmm{12,j} = Mj.peak_x;
        summary_gmm{13,j} = Mj.Kinfl;
        summary_gmm{14,j} = Mj.inf_x;
        summary_gmm{15,j} = Mj.Ipeak;
        summary_gmm{16,j} = Mj.Iinf;
        summary_gmm{17,j} = Mj.g.mu;
        summary_gmm{18,j} = Mj.g.Sigma;
        summary_gmm{19,j} = Mj.weights;
        summary_gmm{20,j} = Mj.Mcompon;
        summary_gmm{21,j} = Mj.polyID;


        medDist{1,j} = Mj.med_dist;
        medDist{2,j} = Mj.comp_medMax;


        g{j}       = Mj.g;
        g_filt{j}  = Mj.g_filt;
        xg{j}      = Mj.xg;
        pg{j}      = Mj.pg;

        if Mj.interest_med
            coeff_out(coeff_interest) = j;
            Idist_win(coeff_interest) = Mj.Idist;
            medPoly(coeff_interest) = Mj.comp_medMax;
            intrst_med(coeff_interest) = Mj.interest_med;
            coeff_interest =1+coeff_interest;
        end

    end
    medDist = medDist';

    best_values = table(coeff_out',Idist_win',medPoly', ...
        intrst_med', ...
        'VariableNames',{'coefficients',...
        'Idist', ...
        'polynomial Idist', ...
        'Idist max >> med'});
    summary_gmm = summary_gmm';
    summary_out = cell2table(summary_gmm,'VariableNames', {'coeff num', 'kj', ...
        'max kj', 'med kj','poly 2 kj(2 max)','Dij','idist', ...
        'med idist','mpdf mu','mpdf std','kpeaks','peak x', ...
        'kinf','inf x','ipeak','iinf','mu gmm','std gmm','weights','M_comp', ...
        'polyID'});
    
    %     %metric for if unimodal or non
    % [max_idist,pIdistMax] = max(Idist);
    % [max_idistMed,pIdist_Med_Max] = max(Idist_med);
    % 
    % [med_idist,pIdistMed_tot] = median(Idist);
    % [med_idistMed,pIdistMedMed] = median(Idist_med);
    % 
    % if max_idist > 2*med_idist
    % 
    % end
    % 
    % if max(Idist_med) > 2*median(Idist_med)
    % 
    % 
    % end
  %  summary_table = [1:numel(coeffs)
end

function M = GMM_basic1D(x,par)
    % One 1-D feature vector x (n x 1) -> metrics + peak/inflection locations
    x = x(:);

    opts  = statset('MaxIter',2000,'TolFun',1e-6,'Display','off');
    K     = 8;
    Mgrid = 100;                 % 100 for paper try 256 if missing def.

    g = fitgmdist(x, K, 'Replicates', 8, ...
        'CovarianceType','diagonal', ...
        'RegularizationValue', 1e-6, ...    
        'Start','plus', 'Options', opts);

    % Grid + pdf


    mu = g.mu(:);                       
    s  = sqrt(squeeze(g.Sigma(:)));      
    a  = g.ComponentProportion(:);    
    polyID = (1:numel(mu))';   % assign original indices from g

    
    

    keep = a > 0.005;   %way to make sure it is not isolated % 2.5% cutoff
    mu = mu(keep);
    s  = s(keep);
    a  = a(keep);
    polyID = polyID(keep);
    Knew = numel(mu);               % actual surviving components

    [a_sort,a_idx] = sort(a,'descend');
    dropThreshold = 4;
    a_drops = (a_sort(1:end-1)-a_sort(2:end))*100;
    a_cutoff = find(a_drops > dropThreshold,1);

    if isempty(a_cutoff) || isnan(a_cutoff)
        % if no cutoff found, take top 3 instead
        topN = min(3, numel(a_sort));
        a_top_idx = a_idx(1:topN);
        a_top = a_sort(1:topN);
    else
        a_top_idx = a_idx(1:min(a_cutoff, numel(a_sort)));
        a_top = a_sort(1:min(a_cutoff, numel(a_sort)));
    end
    mu_sort = mu(a_top_idx);
    s_sort = s(a_top_idx);
    polyID_M = polyID(a_top_idx);
   % a_top = a_sort(1:3);
    M_pdf = 0;
    % combine biggest functions for consolidation metric

    xg = linspace(min(x), max(x), Mgrid).';
    dx = mean(diff(xg));
    scale_a = a_top/sum(a_top);
    for i =1:length(a_top)
        Npdf_i = normpdf(xg,mu_sort(i),s_sort(i)); %scale w a_top(i)?
        M_pdf = M_pdf + scale_a(i)*Npdf_i;     
     end


    if Knew == 0
        mu_new = mean(x);
        s_new  = std(x) + eps;
        a_new  = 1;
        Knew = 1;
    end
    a_new = a / sum(a);
    s_new = reshape(s.^2, 1, 1, Knew);

    g_filt = gmdistribution(mu,s_new,a_new);


    pg = pdf(g_filt, xg);
    pgmax = max(pg);

    % local maxima for peak detection (alternative) but can return
    % prominence, and input minProminence metric. paper uses crit pts of
    % pdf
    
    % promScale = 0.01;
    % s_comp   = sqrt(squeeze(g.Sigma(:)));
    % s_med    = median(s_comp);                     % typical component width
    % minDist  = max(1, round(0.5 * s_med / max(dx,realmin)));
    % 
    % [pks, locs] = findpeaks(pg, 'MinPeakProminence', promScale*max(pg), ...
    %                               'MinPeakDistance',  minDist);
    % isPeak = false(size(pg)); isPeak(locs) = true;

    % derivative to find critical pts for inflection/peaks
    p1  = gradient(pg, dx);
    p2  = gradient(p1, dx);
    
    s1 = p1(1:end-1);
    s2 = p1(2:end);
    crossPosToNeg = (s1 > 0  & s2 <= 0) | (s1 >= 0 & s2 < 0);

    crossLoc = find(crossPosToNeg);  
    % choose the endpoint closer to the zero 
    pick_loc = abs(p1(crossLoc+1)) < abs(p1(crossLoc));
    idx = crossLoc + pick_loc;

    isPeak = false(size(pg));
    validPeak = false(size(idx));
    if ~isempty(idx)
        validPeak = (p2(idx) < 0);             % concave down
        isPeak(idx(validPeak)) = true;
    end

    % inflections look for sign change or 0
    t1 = p2(1:end-1);
    t2 = p2(2:end);
    infl = (t1 >= 0 & t2 <= 0) | (t1 <= 0 & t2 >= 0);

    isInfl = false(size(pg));
    isInfl(2:end) = infl;                      % mark right endpoint

    Ipeak = sum(pg(isPeak)) / pgmax;
    Iinf  = sum(pg(isInfl)) / pgmax;

    [II,JJ] = find(triu(true(Knew),1));

    mu_M = sum(xg .* M_pdf) / sum(M_pdf);  % weighted mean over xg
    std_M = sqrt(sum((xg - mu_M).^2 .* M_pdf) / sum(M_pdf));
    kj = abs(mu - mu_M)/std_M;
    
    if Knew > 1
        Dij = (abs(mu(II) - mu(JJ))).* sqrt(a(II) .* a(JJ)) ...
                ./sqrt(s(II) .* s(JJ));
        dist_per_comp = accumarray([II; JJ], [Dij; Dij], [K 1], @(v){v}, {});
        %avg_dist = cellfun(@mean, dist_per_comp);
        med_dist = cellfun(@median,dist_per_comp);
        %formally idist = median(dij)
        % [~, comp_max] = max(avg_dist);
        % Idist = mean(maxk(avg_dist,2));
        %[Idist, comp_max] = max(avg_dist);
        [Idist,comp_medMax] = max(med_dist);
    else
        Dij = [];  % nothing to compare if only one component remains
        avg_dist = 0; med_dist = 0;
        Idist_val = 0; Idist_med_val = 0;
        %comp_max = 1; 
        comp_medMax = 1;
    end

    if Idist > 1.55*min(med_dist)
        coeff_oInter = 1;
    else
        coeff_oInter = 0;
    end

    % if Idist_med > 1.55*min(med_dist)
    %     coeff_oInterMed = 1;
    % else
    %     coeff_oInterMed = 0;
    % end
    % Package extras
    peak_x  = xg(isPeak);
    inf_x   = xg(isInfl);
    Kpeaks  = numel(peak_x);
    Kinfl = numel(inf_x);

    %output idist max and poly from max idist

    M = struct('Ipeak',Ipeak, ...
           'Iinf',Iinf, ...
           'Dij',Dij, ...
           'Idist',Idist, ...
           'g',g,'g_filt',g_filt, 'xg',xg, 'pg',pg, ...
           'Kpeaks',Kpeaks, 'peak_x',peak_x, ...
           'inf_x',inf_x,'Kinfl',Kinfl, ...
            'med_dist',med_dist, ...
            'comp_medMax',comp_medMax, ...
            'interest_med',coeff_oInter, ...
            'dist_kj',kj,'mu_mpdf',mu_M,'s_mpdf',std_M,'weights',a_new,'Mcompon', polyID_M, ...
            'polyID',polyID);
end


function Xwav = haar_coeffs(spikes,par)
    scales = par.scales;
    feature = par.features;
    nspk = size(spikes,1);
    ls = size(spikes,2);
    cc=zeros(nspk,ls);
		try
			spikes_l = reshape(spikes',numel(spikes),1);
			if exist('wavedec')
				[c_l,l_wc] = wavedec(spikes_l,scales,'haar');
			else
				[c_l,l_wc] = fix_wavedec(spikes_l,scales);
			end
			wv_c = [0;l_wc(1:end-1)];
			nc = wv_c/nspk;
			wccum = cumsum(wv_c);
			nccum = cumsum(nc);
			for cf = 2:length(nc)
				cc(:,nccum(cf-1)+1:nccum(cf)) = reshape(c_l(wccum(cf-1)+1:wccum(cf)),nc(cf),nspk)';
			end
		catch
		    if exist('wavedec')                             % Looks for Wavelets Toolbox
				for i=1:nspk                                % Wavelet decomposition
					[c,l] = wavedec(spikes(i,:),scales,'haar');
					cc(i,1:ls) = c(1:ls);
				end
			else
				for i=1:nspk                                % Replaces Wavelets Toolbox, if not available
					[c,l] = fix_wavedec(spikes(i,:),scales);
					cc(i,1:ls) = c(1:ls);
				end
			end
		
		end
		
        % mask = true(size(cc,1),1);
        % for i = 1:ls
        %     mu  = mean(cc(:,i));
        %     sd  = std(cc(:,i));
        %     thr = 3*sd;
        %     lo  = mu - thr;
        %     hi  = mu + thr;
        %     mask = mask & (cc(:,i) >= lo) & (cc(:,i) <= hi);
        % end
        %Xwav = cc(mask, :);
        Xwav = cc;
end

function [mu_new,s_new,a_new,polyID_new] = mergeNearby(mu,s,a,polyID,tol)
    % pairwise distance between means
    n = numel(mu);
    merged = false(1,n);
    mu_new = []; s_new = []; a_new = []; polyID_new = {};

    for i = 1:n
        if merged(i), continue; end
        % find group close to mu(i)
        group = find(abs(mu - mu(i)) < tol*mean(s));
        merged(group) = true;

        % combine weights
        a_sum = sum(a(group));
        % weighted mean
        mu_c = sum(a(group).*mu(group)) / a_sum;
        % pooled variance
        var_c = sum(a(group).*(s(group).^2 + mu(group).^2))/a_sum - mu_c^2;
        s_c = sqrt(max(var_c,eps));

        % store
        mu_new(end+1,1) = mu_c;
        s_new(end+1,1)  = s_c;
        a_new(end+1,1)  = a_sum;
        polyID_new{end+1,1} = polyID(group);
    end
end


