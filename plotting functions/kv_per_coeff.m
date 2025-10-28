function kv_per_coeff(ks_out,folderName,channelNum,maxK_sort)
    % change for idist over ks over lay ii double plot
    % need to lable
    fig_k = figure;
    plot(maxK_sort(:,2))
    title(" max k vals per coeff")
    for k = 1:length(maxK_sort(:,2))
       text(k+0.1,maxK_sort(k,2)+0.1,num2str(maxK_sort(k,1))); 
    end
    filename_kv = fullfile(folderName,sprintf('ch%s_maxKvalsPerCoef.png', channelNum));
    exportgraphics(fig_k, filename_kv, 'Resolution', 300);

end