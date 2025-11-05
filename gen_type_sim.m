function gen_type_sim
    % this function takes in real data and makes it in the form of 
    % of sim data - mat file w spikes and cluster_class
    [files, path] = uigetfile('*.mat', 'Select .mat files', 'MultiSelect', 'on');
    if ischar(files)
        files = {files};
    end

    % match by test/channel num
    segMentFile = split(files,["_","."]);
    segMentFile = squeeze(segMentFile);
    % see if times file or not
    times_files_idx = find(segMentFile(:,1) == "times");
    spikes_files_idx = find(segMentFile(:,1) ~= "times");


    filesTimes = files(:,times_files_idx);

    filesSpikes = files(:,spikes_files_idx);

    segMent_times = segMentFile(times_files_idx,:);

    spike_segment = segMentFile(spikes_files_idx,:);

    for i = 1:length(filesTimes)
        fullTimes = fullfile(path, filesTimes{i});  
        
        % find equivalent spike (matching test and channel)
        match_idx = ismember(spike_segment(:,2:3), segMent_times(i,2:3), 'rows');
        if ~any(match_idx)
            warning('No matching spike file for %s', filesTimes{i});
            continue;
        end
        
        spike_full = fullfile(path, filesSpikes{find(match_idx, 1)}); 
    
        % Load data
        times = load(fullTimes);
        test_data = load(spike_full);
    
        % Save new combined file
        cluster_class = times.cluster_class;
        spikes = test_data.spikes;
    
        dataName = sprintf('sim_%s_%s.mat', segMent_times{i,2}, segMent_times{i,3});
        save(fullfile(path, dataName), 'spikes', 'cluster_class');
    end

end