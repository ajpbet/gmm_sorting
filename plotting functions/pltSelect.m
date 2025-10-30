function pltSelect(len_select,vector,channelNum,folderName,plotID)
    % input variables
    % select is var selected by knee [value,coeff#,gauss#]
    % len_select is how long select vector is
    % vector is full sorted vector [value,coeff#,gauss#]

     % Create figure
    fig = figure('Visible', 'on'); % Change to 'off' for batch mode
    hold on
    grid on
    
    % plot type
    switch plotID
        case 1
            plot_type = 'kv';
        case 2
            plot_type = 'kvNoPkExc';
        case 3
            plot_type = 'medDist';
        case 4
            plot_type = 'medDistNoPkExc';
    end
        
    title(sprintf('%s per Coefficient',plot_type));
    xlabel('Coefficient Index (Sorted)');
    ylabel(sprintf('%s Value',plot_type));
    plot(vector(:,1), 'b-', 'DisplayName', plot_type);
    
    plot_val = length(vector(:,1)) -len_select;
    % Plot IDist points and markers
    for k = 1:length(vector(:,1))
        if k >= plot_val
            plot(k, vector(k,1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
        else
            plot(k, vector(k,1), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end

    % Label coefficient numbers for IDist
    for k = 1:length(vector(:,1))
        text(k + 0.1, vector(k,1) + 0.1, num2str(vector(k,2)), 'Color', 'b', ...
             'FontSize', 8, 'HorizontalAlignment', 'left');
        text(k + 0.1, vector(k,1) + 0.15, sprintf('G:%d',vector(k,3)), 'Color', 'b', ...
             'FontSize', 8, 'HorizontalAlignment', 'left');
    end


    % Add reference lines
    xline(plot_val, '--', 'Color', 'red', 'LineWidth', 0.5, 'Alpha', 0.7, 'DisplayName', 'knee cutoff');

    % Build legend
    h_match(1) = plot(NaN, NaN, 'b-', 'DisplayName', plot_type);
    h_match(2) = scatter(NaN, NaN, 80, 'k', 'filled', 'DisplayName', 'Selected');
    h_match(3) = scatter(NaN, NaN, 80, 'kx', 'DisplayName', 'Not selected');
    h_match(4) = plot(NaN, NaN, 'r-', 'DisplayName', 'knee cutoff');

    legend(h_match, 'Location', 'best');
    hold off;

    % Format output filename and save
    if isnumeric(channelNum)
        channelNum = num2str(channelNum);
    end

    filename = fullfile(folderName, sprintf('ch%s_plt_%s.png', channelNum,plot_type));
    exportgraphics(fig, filename, 'Resolution', 300);

end