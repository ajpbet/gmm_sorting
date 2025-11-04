function sim_run
    % file selection
    [files, path] = uigetfile('*.mat', 'Select .mat files', 'MultiSelect', 'on');
    if ischar(files)
        files = {files};
    end
    fullPaths = fullfile(path, files);
    
    num_units = zeros(1, numel(fullPaths));

    % Folder setup
    folderName = 'results';
    folderGMM  = fullfile(folderName, 'GMM');
    folderCoeff = fullfile(folderName, 'Coeff');
    if ~exist(folderName, 'dir'), mkdir(folderName); end
    if ~exist(folderGMM, 'dir'), mkdir(folderGMM); end
    if ~exist(folderCoeff, 'dir'), mkdir(folderCoeff); end

    % === Summary structure (NEW) ===
    coeffSummary = struct('file', [], 'select_all', [], 'ks_coeff', [], ...
        'select_spike_match', []);

    par = set_parameters;
    summaryFile = fullfile(folderName, 'coeff_summary_all.mat');
    % loop through files
    for i = 1:numel(fullPaths)
        data = load(fullPaths{i});
        
        spikes = data.spikes;
        cluster_class = data.cluster_class;
        num_units(i) = numel(unique(cluster_class)) - 1;

        % Base name
        [~, baseName, ~] = fileparts(fullPaths{i});

        % --- GMM File Path ---
        gmmFileName = fullfile(folderGMM, ['gmm_' baseName '.mat']);

        % --- Run GMM Extract ---
        [select_all, ks_coeff,select_spike_match] = ...
            GMM_extract(spikes, cluster_class, par, gmmFileName,baseName);

        % --- Coeff File Path ---
        coeffFileName = fullfile(folderCoeff, ['coeff_' baseName '.mat']);
        save(coeffFileName, 'select_all', 'ks_coeff','select_spike_match');

        % --- Add to Summary (NEW) ---
        coeffSummary(i).file = baseName;
        coeffSummary(i).select_all = select_all;
      %  coeffSummary(i).select_allNoPk = select_allNoPk;
        coeffSummary(i).ks_coeff = ks_coeff;
        coeffSummary(i).select_spike_match = select_spike_match;

        fprintf('Saved:\n  %s\n  %s\n', gmmFileName, coeffFileName);
    end
    save(summaryFile, 'coeffSummary');
    fprintf('Saved summary: %s\n', summaryFile);
% else
%     for i = 1:numel(fullPaths)
%         data = load(fullPaths{i});
%         cluster_class = data.cluster_class;
%         num_units(i) = numel(unique(cluster_class)) - 1;
%     end
%         load(summaryFile);
    %end

    %% plotting functions
    nFiles = numel(coeffSummary);  % 95
    
    % Preallocate cell arrays
    select_gauss0    = cell(nFiles,1);
    select_gauss1pct = cell(nFiles,1);
    select_gauss2_5pct = cell(nFiles,1);
    select_gauss5pct = cell(nFiles,1);
    select_gauss10pct = cell(nFiles,1);
    ks_coeff_all    = cell(nFiles,1);
    
    select_spikeMch = cell(nFiles,1);

    for i = 1:nFiles
        % Each unit's select_all field contains the 5 structs
        sel_all = coeffSummary(i).select_all;
        ks_coeff_all{i} = coeffSummary(i).ks_coeff;
        select_gauss0{i}     = sel_all.original.select_gauss(:,3);
        select_gauss1pct{i}  = sel_all.trim1_0pct.select_gauss(:,3);
        select_gauss2_5pct{i} = sel_all.trim2_5pct.select_gauss(:,3);
        select_gauss5pct{i}  = sel_all.trim5_0pct.select_gauss(:,3);
        select_gauss10pct{i} = sel_all.trim10_0pct.select_gauss(:,3);

        select_spikeMch{i} = coeffSummary(i).select_spike_match(:,3);
    end

    % plot ks vs gmm coeff matches
    plotKSCoeffOverlap(ks_coeff_all, select_gauss0, num_units, 'overlap ks v. gmm 0% cutoff',1);
    plotKSCoeffOverlap(ks_coeff_all, select_gauss1pct, num_units, 'overlap ks v. gmm 1% cutoff',2);
    plotKSCoeffOverlap(ks_coeff_all, select_gauss2_5pct, num_units, 'overlap ks v. gmm 2.5% cutoff',3);
    plotKSCoeffOverlap(ks_coeff_all, select_gauss5pct, num_units, 'overlap ks v. gmm 5% cutoff',4);
    plotKSCoeffOverlap(ks_coeff_all, select_gauss10pct, num_units, 'overlap ks v. gmm 10% cutoff',5);
    plotKSCoeffOverlap(ks_coeff_all, select_spikeMch, num_units, 'overlap ks v. spike match 0%',5);

    % plot ks vs gmm coeff mismatches

    plotKSCoeffMismatch(ks_coeff_all, select_gauss0, num_units, 'mismatch ks v. gmm 0% cutoff',1);
    plotKSCoeffMismatch(ks_coeff_all, select_gauss1pct, num_units, 'mismatch ks v. gmm 1% cutoff',2);
    plotKSCoeffMismatch(ks_coeff_all, select_gauss2_5pct, num_units, 'mismatch ks v. gmm 2.5% cutoff',3);
    plotKSCoeffMismatch(ks_coeff_all, select_gauss5pct, num_units, 'mismatch ks v. gmm 5% cutoff',4);
    plotKSCoeffMismatch(ks_coeff_all, select_gauss10pct, num_units, 'mismatch ks v. gmm 10% cutoff',5);
    plotKSCoeffMismatch(ks_coeff_all, select_spikeMch, num_units, 'mismatch ks v. spike match',5);


    % plot ks and gmm  coeff nums vs unit numbers 
    plotKSCoeffCount(ks_coeff_all, select_gauss0, num_units, 'coeff ct ks v. gmm 0% cutoff',1);
    plotKSCoeffCount(ks_coeff_all, select_gauss1pct, num_units, 'coeff ct ks v. gmm 1% cutoff',2);
    plotKSCoeffCount(ks_coeff_all, select_gauss2_5pct, num_units, 'coeff ct ks v. gmm 2.5% cutoff',3);
    plotKSCoeffCount(ks_coeff_all, select_gauss5pct, num_units, 'coeff ct ks v. gmm 5% cutoff',4);
    plotKSCoeffCount(ks_coeff_all, select_gauss10pct, num_units, 'coeff ct ks v. gmm 10% cutoff',5);
    plotKSCoeffCount(ks_coeff_all, select_spikeMch, num_units, 'coeff ct ks v. gmm spike match',5);

end
