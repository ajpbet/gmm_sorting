function [top5_iPk,top5_iDist,top5_iInf,mse_vals,table_summary_detect,table_summary_vals,table_combined] = plot_results(pg,xg,wd_coeff,ks_coeff,iPk,iDist,iInf,g,kpeaks,kinfl)
    [desc_iPk,orgIpkloc] = sort(iPk,'descend');
    top5_iPk = orgIpkloc(1:5);
    
    [desc_iDist,orgIDloc] = sort(iDist,'descend');
    top5_iDist = orgIDloc(1:5);

    [desc_iInf,orgIinfloc] = sort(iInf,'descend');
    top5_iInf = orgIinfloc(1:5);

    all_coeff = [top5_iPk;top5_iDist;top5_iInf;ks_coeff'];

    foundWD = unique(all_coeff);

    iPk_present = zeros(length(foundWD),1);
    iDist_present = zeros(length(foundWD),1);
    iInf_present = zeros(length(foundWD),1);
    ks_present = zeros(length(foundWD),1);

    [kde_pdf,kde_xf] = kde_est(wd_coeff,foundWD);

    mse_vals = mean_square_error(pg,xg,kde_pdf,kde_xf,foundWD);

    for i = 1:length(foundWD)
        iPk_present(i) = ismember(foundWD(i),top5_iPk);
        iDist_present(i) = ismember(foundWD(i),top5_iDist);
        iInf_present(i) = ismember(foundWD(i),top5_iInf);
        ks_present(i) = ismember(foundWD(i),ks_coeff);


        if iPk_present(i) 
            iPk_title = " *Ipk* = " + iPk(foundWD(i)) + ", ";
        else
            iPk_title = "Ipk = " + iPk(foundWD(i)) + ", ";
        end
        if iDist_present(i) 
            iDist_title = "*IDist* = " + iDist(foundWD(i)) + ", ";
        else
            iDist_title = "IDist = " + iDist(foundWD(i)) + ", ";
        end

        if iInf_present(i) 
            iInf_title = "*Iinf* = " + iInf(foundWD(i)) + ", ";
        else
            iInf_title = "Iinf = " + iInf(foundWD(i)) + ", ";
        end

        if ks_present(i) 
            ks_title = "*KS* ";
        else
            ks_title = " ";
        end

        figure;
        subplot(2,1,1)
        histogram(wd_coeff(:,foundWD(i)),1000,'Normalization','pdf');
        title("coeff " + foundWD(i) + " for: " + iPk_title + iDist_title + iInf_title + ks_title);
        hold on
        plot(kde_xf{:,i},kde_pdf{:,i},'LineWidth',2);
        legend('hist','kde')
        hold off

        subplot(2,1,2)
        plot(xg{:,foundWD(i)},pg{:,foundWD(i)},'k-','LineWidth',2); hold on

        gmm_model = g{foundWD(i)};
        mu = gmm_model.mu(:);
        sigma = sqrt(squeeze(gmm_model.Sigma(:)));
        w = gmm_model.ComponentProportion(:);
        cmap = lines(length(mu));

        for k = 1:length(mu)
            comp_pdf = gmm_model.ComponentProportion(k) * normpdf(xg{:,foundWD(i)}, mu(k), sigma(k));
            gmm_model.ComponentProportion(k);
            plot(xg{:,foundWD(i)}, comp_pdf,'Color',[0.7 0.7 0.7],'LineWidth',1);
            [~,pk_idx] = max(comp_pdf);
            x_vals = xg{foundWD(i)};
            plot(x_vals(pk_idx), comp_pdf(pk_idx), 'o', 'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 6);
            text(x_vals(pk_idx), comp_pdf(pk_idx), sprintf('plot=%.2f', k), 'Color', cmap(k,:), 'VerticalAlignment', 'bottom',...
                'HorizontalAlignment', 'center');

        end
        text(mean(xg{:,foundWD(i)}), max(pg{:,foundWD(i)})*1.05, sprintf('MSE = %.4f', mse_vals(i)),...
            'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');

        title("pdf for coeff: " + foundWD(i));

        plot(kde_xf{:,i},0.01 + kde_pdf{:,i},'LineWidth',2);
        legend('pdf','','','','','',...
            '','','','','','','','',...
            '','','','kde');
        hold off

        xlims = xlim; ylims = ylim;
        info_str = cell(length(mu)+1,1);
        info_str{1} = sprintf('Comp | Mean | Std | Weight');
        for k = 1:length(mu)
            info_str{k+1} = sprintf('%3d | %.2f | %.2f | %.2f', k, mu(k), sigma(k), w(k));
        end

        text(xlims(2)*0.6, mean(ylims), info_str,'Units','data','HorizontalAlignment','left','FontName','Courier','FontSize',9);
    end
    table_summary_detect = table(foundWD, mse_vals', iPk_present, iDist_present, iInf_present, ks_present, ...
        (iPk_present+iDist_present+iInf_present+ks_present),...
        'VariableNames', {'coeffs','mseVals','IPeak_Detect','IDist_Detect','IInfl_Detect',...
        'KS_Detect','sum detect'});

    iPeak_rel = iPk(foundWD);
    iDist_rel = iDist(foundWD);
    iInf_rel = iInf(foundWD);
    kPeak_rel = kpeaks(foundWD);
    kInfl_rel = kinfl(foundWD);

    table_summary_vals = table(foundWD, mse_vals', iPeak_rel, iDist_rel, iInf_rel,...
        kPeak_rel, kInfl_rel,...
        'VariableNames', {'coeffs','mseVals','iPeak','iDist','iInf',...
        'Kpeaks','Kinfl'});

    iPeak_rel(~iPk_present) = 0;
    iDist_rel(~iDist_present) = 0;
    iInf_rel(~iInf_present) = 0;


    table_combined = table(foundWD, mse_vals', iPeak_rel, iDist_rel, iInf_rel,ks_present, ...
        kPeak_rel, kInfl_rel,...
        'VariableNames', {'coeffs','mseVals','iPeak','iDist','iInf','KS Present','Kpeaks',...
        'Kinfl'});

end


function [kde_pdf,kde_xf] = kde_est(wd_coeff,foundWD)
    kde_pdf = cell(1, length(foundWD));
    kde_xf  = cell(1, length(foundWD));
    for i = 1:length(foundWD)
        [kde_pdf{i},kde_xf{i}] = kde(wd_coeff(:,foundWD(i)));
    end
end

% 
function mse_vals = mean_square_error(pg,xg,kde_pdf,kde_xf,foundWD)
    for i = 1:length(foundWD)
        kde_pdf_interp = interp1(kde_xf{i}, kde_pdf{i}, xg{foundWD(i)}, 'linear');
        mse_vals(i) = mean((normalize(kde_pdf_interp) - normalize(pg{foundWD(i)})).^2);        
    end

end