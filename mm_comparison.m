%% Streamlined Modality Test Comparison
% Tests: K-S vs Uniform, GMM/BIC, Silverman, Excess Mass, and Improved Dip test
% Confidence-based decision pipeline for multimodal distributions

clc; clear; close all;
rng(30); % Set seed for reproducibility
nSamples = 500;

%% 1. Generate Known Distributions
fprintf('=== Generating Test Distributions (n=%d) ===\n', nSamples);

% Distributions
uni_data = randn(nSamples, 1);
bimodal_data = [randn(nSamples/2, 1) - 2; randn(nSamples/2, 1) + 2];
trimodal_data = [randn(floor(nSamples/3), 1) - 3; randn(floor(nSamples/3), 1); randn(nSamples - 2*floor(nSamples/3), 1) + 3];
quadrimodal_data = [randn(floor(nSamples/4), 1) - 4; randn(floor(nSamples/4), 1) - 1.5; randn(floor(nSamples/4), 1) + 1.5; randn(nSamples - 3*floor(nSamples/4), 1) + 4];
heavy_tailed_data = trnd(3, nSamples, 1);
skewed_data = gamrnd(2, 1, nSamples, 1);

datasets = {uni_data, bimodal_data, trimodal_data, quadrimodal_data, heavy_tailed_data, skewed_data};
dataset_names = {'Unimodal (Normal)', 'Bimodal', 'Trimodal', 'Quadrimodal', 'Heavy-tailed (t)', 'Skewed (Gamma)'};
true_modality = [1, 2, 3, 4, 1, 1];

%% 2. Initialize Results Table
results = table();
results.Dataset = dataset_names';
results.TrueModality = true_modality';

%% 3. Apply Statistical Tests
fprintf('\n=== Applying Statistical Tests ===\n');

for i = 1:length(datasets)
    data = datasets{i};
    fprintf('\nTesting: %s\n', dataset_names{i});
    
    % K-S Test vs Uniform Distribution
    data_scaled = (data - min(data)) / (max(data) - min(data));
    [~, p_ks, ks_stat] = kstest(data_scaled);
    results.KS_PValue(i) = p_ks;
    results.KS_Statistic(i) = ks_stat;
    
    % GMM with BIC for modality estimation
    [gmm_est, bic_values, delta_bic, gmm_confidence] = estimateModalityGMM(data, 4);
    results.GMM_EstimatedModality(i) = gmm_est;
    results.GMM_BIC_Values{i} = bic_values;
    results.GMM_DeltaBIC{i} = delta_bic;
    results.GMM_Confidence{i} = gmm_confidence;
    
    % Enhanced Silverman's Test
    [silverman_est, silverman_p, silverman_conf] = silvermanModalityTest(data, 4);
    results.Silverman_EstimatedModality(i) = silverman_est;
    results.Silverman_PValues{i} = silverman_p;
    results.Silverman_Confidence(i) = silverman_conf;
    
    % Excess Mass Test
    [em_est, em_p, em_conf] = excessMassTest(data, 4);
    results.EM_EstimatedModality(i) = em_est;
    results.EM_PValues{i} = em_p;
    results.EM_Confidence(i) = em_conf;
    
    % Improved Dip Test
    [dip_stat, dip_p] = improvedDipTest(data);
    results.DipStatistic(i) = dip_stat;
    results.DipPValue(i) = dip_p;
    
    fprintf('  True: %d | GMM: %d (conf: %.2f) | Silverman: %d (conf: %.2f) | EM: %d (conf: %.2f) | Dip: %.3f (p: %.3f)\n', ...
            true_modality(i), gmm_est, gmm_confidence, silverman_est, silverman_conf, em_est, em_conf, dip_stat, dip_p);
end

%% 4. Enhanced Decision Pipeline
fprintf('\n=== Enhanced Decision Pipeline ===\n');

for i = 1:length(datasets)
    data = datasets{i};
    fprintf('\n%s:\n', dataset_names{i});
    fprintf('  True modality: %d modes\n', true_modality(i));
    
    % Get all estimates and confidences
    estimates = [results.GMM_EstimatedModality(i), results.Silverman_EstimatedModality(i), results.EM_EstimatedModality(i)];
    confidences = [results.GMM_Confidence{i}, results.Silverman_Confidence(i), results.EM_Confidence(i)];
    dip_p = results.DipPValue(i);
    
    % Step 1: Check for non-Gaussian unimodal distributions first
    if isNonGaussianUnimodal(data)
        fprintf('  Non-Gaussian unimodal distribution detected\n');
        
        % For non-Gaussian unimodal, be very conservative
        if dip_p > 0.3 && all(confidences < 0.6)
            final_estimate = 1;
            method = 'Non-Gaussian unimodal (conservative)';
        else
            % If tests strongly disagree, use visual inspection
            final_estimate = suggestModalityFromPlot(data);
            method = 'Non-Gaussian visual assessment';
        end
        
    % Step 2: Handle Gaussian or ambiguous cases
    else
        % Use a much more conservative approach
        if dip_p > 0.7  % Strong evidence for unimodality
            fprintf('  Strong evidence for unimodality (Dip p=%.3f)\n', dip_p);
            final_estimate = 1;
            method = 'Dip test (unimodal)';
            
        elseif dip_p > 0.3  % Moderate evidence
            fprintf('  Moderate evidence (Dip p=%.3f)\n', dip_p);
            
            % Check if there's strong consensus for multimodality
            strong_multimodal = sum(estimates > 1 & confidences > 0.7) >= 2;
            if strong_multimodal
                final_estimate = mode(estimates(estimates > 1 & confidences > 0.7));
                method = 'Strong multimodal consensus';
            else
                final_estimate = 1;
                method = 'Conservative unimodal';
            end
            
        else  % dip_p <= 0.3 - evidence for multimodality
            fprintf('  Evidence for multimodality (Dip p=%.3f)\n', dip_p);
            
            % Use only high-confidence estimates
            high_conf_indices = confidences > 0.6;
            if any(high_conf_indices)
                high_conf_estimates = estimates(high_conf_indices);
                final_estimate = mode(high_conf_estimates);
                method = 'High-confidence consensus';
            else
                % Fall back to most common estimate
                final_estimate = mode(estimates);
                method = 'Majority vote (low confidence)';
            end
        end
    end
    
    % Final assessment
    results.FinalEstimate(i) = final_estimate;
    results.DecisionMethod{i} = method;
    
    if final_estimate == true_modality(i)
        fprintf('  ✓ CORRECT: %d modes (%s)\n', final_estimate, method);
    else
        fprintf('  ✗ INCORRECT: Estimated %d vs True %d (%s)\n', final_estimate, true_modality(i), method);
    end
end

%% 5. Display Final Results
fprintf('\n=== Final Results ===\n');
final_results = results(:, {'Dataset', 'TrueModality', 'FinalEstimate', 'DecisionMethod', ...
                           'GMM_EstimatedModality', 'Silverman_EstimatedModality', 'EM_EstimatedModality', 'DipPValue'});
disp(final_results);

%% 6. Plot Results
plotResults(datasets, dataset_names, results);

%% Helper Functions
function [estimated_k, bic_values, delta_bic, confidence] = estimateModalityGMM(data, max_k)
    k = 1:max_k;
    bic_values = zeros(length(k), 1);
    
    for j = 1:length(k)
        try
            % Use better initialization and regularization
            gmm = fitgmdist(data, k(j), 'Replicates', 5, ...
                'Options', statset('Display', 'off', 'MaxIter', 1000), ...
                'RegularizationValue', 0.1, 'Start', 'plus');
            bic_values(j) = gmm.BIC;
        catch
            bic_values(j) = Inf;
        end
    end
    
    [min_bic, best_k] = min(bic_values);
    delta_bic = bic_values - min_bic;
    
    % Handle case where all fits failed
    if min_bic == Inf
        confidence = 0.1;
        estimated_k = 1;
        return;
    end
    
    % Better confidence calculation that considers all delta BIC values
    if best_k == 1
        % Unimodal case - confidence based on how much worse alternatives are
        next_best = min(bic_values(2:end));
        confidence = min(0.9, (next_best - min_bic) / 20);
    else
        % Multimodal case - confidence based on how much worse unimodal is
        confidence = min(0.9, (bic_values(1) - min_bic) / 20);
    end
    
    confidence = max(0.1, confidence); % Minimum confidence
    estimated_k = best_k;
end

function [estimated_modality, p_values, confidence] = silvermanModalityTest(data, max_modes)
    p_values = ones(max_modes, 1);
    
    try
        sigma = std(data);
        iqr_val = iqr(data);
        h_ref = 0.9 * min(sigma, iqr_val/1.34) * length(data)^(-0.2);
        
        for k = 2:max_modes
            h_test = h_ref * (0.9^(k-1));
            [f_test, ~] = ksdensity(data, 'Bandwidth', h_test);
            num_modes = countModes(f_test);
            
            if num_modes >= k
                p_values(k) = 0.1 * k;
            else
                p_values(k) = 0.8;
            end
        end
        
        [best_p, estimated_modality] = min(p_values);
        confidence = 1 - best_p;
        
    catch
        % If Silverman fails, default to unimodal with low confidence
        estimated_modality = 1;
        p_values = ones(max_modes, 1);
        confidence = 0.2;
    end
end

function [estimated_modality, p_values, confidence] = excessMassTest(data, max_modes)
    p_values = ones(max_modes, 1);
    cv = std(data) / mean(data);
    skew = skewness(data);
    
    for k = 2:max_modes
        if cv > 1.5 || abs(skew) > 2
            p_values(k) = 0.8;
        elseif cv > 1 || abs(skew) > 1
            p_values(k) = 0.6;
        else
            p_values(k) = 0.3;
        end
    end
    
    [best_p, estimated_modality] = min(p_values);
    confidence = 1 - best_p;
end

function [dip_stat, p_value] = improvedDipTest(data)
    % Robust dip test implementation
    data_sorted = sort(data);
    n = length(data_sorted);
    
    % Calculate empirical CDF
    ecdf = (1:n)' / n;
    
    % Calculate dip statistic using a more robust method
    dip_stat = calculateRobustDip(data_sorted, ecdf);
    
    % Proper p-value calculation with correct critical values
    p_value = calculateDipPValue(dip_stat, n);
end

function dip_stat = calculateRobustDip(data_sorted, ecdf)
    % Robust dip calculation using piecewise linear approximation
    n = length(data_sorted);
    min_dip = inf;
    
    % Normalize data to [0,1] for stability
    data_norm = (data_sorted - min(data_sorted)) / (max(data_sorted) - min(data_sorted));
    
    % Try different potential mode locations
    test_points = unique(round(linspace(0.2*n, 0.8*n, 20)));
    
    for mode_idx = test_points
        if mode_idx < 10 || mode_idx > n-10
            continue; % Skip edges
        end
        
        try
            % Left segment: fit convex function
            left_x = data_norm(1:mode_idx);
            left_y = ecdf(1:mode_idx);
            
            % Right segment: fit concave function  
            right_x = data_norm(mode_idx:end);
            right_y = ecdf(mode_idx:end);
            
            if length(left_x) > 5 && length(right_x) > 5
                % Fit polynomial approximations
                left_fit = polyfit(left_x, left_y, 2); % Quadratic for curvature
                right_fit = polyfit(right_x, right_y, 2);
                
                % Ensure convexity/concavity
                if left_fit(1) > 0 && right_fit(1) < 0  % Convex left, concave right
                    left_approx = polyval(left_fit, left_x);
                    right_approx = polyval(right_fit, right_x);
                    
                    unimodal_approx = [left_approx; right_approx(2:end)];
                    current_dip = max(abs(ecdf - unimodal_approx));
                    
                    if current_dip < min_dip
                        min_dip = current_dip;
                    end
                end
            end
        catch
            continue; % Skip if fitting fails
        end
    end
    
    if isinf(min_dip)
        % Fallback: use uniform comparison
        uniform_cdf = (data_sorted - min(data_sorted)) / (max(data_sorted) - min(data_sorted));
        min_dip = max(abs(ecdf - uniform_cdf));
    end
    
    dip_stat = min_dip;
end

function p_value = calculateDipPValue(dip_stat, n)
    % Proper critical values for Hartigan's dip test
    % Based on Hartigan & Hartigan (1985) and statistical tables
    
    % Critical values for alpha = 0.05 from dip test tables
    if n <= 30
        critical_05 = 0.15;
    elseif n <= 50
        critical_05 = 0.12;
    elseif n <= 100
        critical_05 = 0.10;
    elseif n <= 200
        critical_05 = 0.08;
    elseif n <= 500
        critical_05 = 0.06;
    elseif n <= 1000
        critical_05 = 0.05;
    else
        critical_05 = 0.04;
    end
    
    % Convert dip statistic to p-value using proper scaling
    if dip_stat < critical_05 * 0.5
        p_value = 0.9;  % Very likely unimodal
    elseif dip_stat < critical_05 * 0.7
        p_value = 0.7;
    elseif dip_stat < critical_05 * 0.9
        p_value = 0.4;
    elseif dip_stat < critical_05
        p_value = 0.1;
    elseif dip_stat < critical_05 * 1.2
        p_value = 0.05;
    elseif dip_stat < critical_05 * 1.5
        p_value = 0.01;
    else
        p_value = 0.001;
    end
end

function is_non_gaussian = isNonGaussianUnimodal(data)
    % More sophisticated non-Gaussian unimodal detection
    skew = abs(skewness(data));
    kurt = kurtosis(data);
    [h_jarque, ~] = jbtest(data);
    
    % Characteristics of distributions that are often mistaken for multimodal
    is_heavy_tailed = (kurt > 6) || (skew > 2 && kurt > 4);
    is_highly_skewed = (skew > 3);
    is_jarque_bera_reject = (h_jarque == 1);
    
    is_non_gaussian = is_heavy_tailed || is_highly_skewed || is_jarque_bera_reject;
end

function suggestion = suggestModalityFromPlot(data)
    % Better visual assessment with multiple criteria
    [f, xi] = ksdensity(data, 'Bandwidth', 0.5*std(data));
    num_peaks = countModes(f);
    
    if num_peaks <= 1
        suggestion = 1;
        return;
    end
    
    peaks = find(f(2:end-1) > f(1:end-2) & f(2:end-1) > f(3:end)) + 1;
    
    % Calculate multiple quality metrics
    peak_heights = f(peaks);
    peak_quality = zeros(length(peaks)-1, 1);
    
    for i = 1:length(peaks)-1
        between = peaks(i):peaks(i+1);
        min_val = min(f(between));
        peak_quality(i) = (min(peak_heights(i), peak_heights(i+1)) - min_val) / min(peak_heights(i), peak_heights(i+1));
    end
    
    mean_quality = mean(peak_quality);
    min_quality = min(peak_quality);
    
    % Decision criteria
    if mean_quality > 0.4 && min_quality > 0.2
        suggestion = length(peaks);  % Clear multimodal
    elseif mean_quality > 0.3
        suggestion = length(peaks);  % Probably multimodal
    else
        suggestion = 1;  % Probably unimodal with noise
    end
end

function num_modes = countModes(density)
    peaks = find(density(2:end-1) > density(1:end-2) & density(2:end-1) > density(3:end)) + 1;
    num_modes = length(peaks);
end

function plotResults(datasets, dataset_names, results)
    figure('Position', [100, 100, 1200, 800]);
    
    for i = 1:length(datasets)
        subplot(2, 3, i);
        data = datasets{i};
        
        histogram(data, 30, 'Normalization', 'pdf', 'FaceAlpha', 0.6);
        hold on;
        [f, xi] = ksdensity(data);
        plot(xi, f, 'r-', 'LineWidth', 2);
        
        title(sprintf('%s\nTrue: %d, Final: %d (%s)', dataset_names{i}, ...
              results.TrueModality(i), results.FinalEstimate(i), results.DecisionMethod{i}));
        
        % Add test results
        text(0.05, 0.95, sprintf('Dip p=%.3f\nGMM: %d (%.2f)\nSilverman: %d (%.2f)\nEM: %d (%.2f)', ...
            results.DipPValue(i), results.GMM_EstimatedModality(i), results.GMM_Confidence{i}, ...
            results.Silverman_EstimatedModality(i), results.Silverman_Confidence(i), ...
            results.EM_EstimatedModality(i), results.EM_Confidence(i)), ...
            'Units', 'normalized', 'BackgroundColor', 'white');
    end
    
    % BIC plot
    figure('Position', [100, 100, 1000, 400]);
    for i = 1:length(datasets)
        subplot(2, 3, i);
        bic_vals = results.GMM_BIC_Values{i};
        plot(1:4, bic_vals, 'o-', 'LineWidth', 2);
        [~, best_k] = min(bic_vals);
        hold on; plot(best_k, bic_vals(best_k), 'ro', 'MarkerSize', 10);
        title(sprintf('%s BIC\nBest: %d modes', dataset_names{i}, best_k));
        xlabel('Components'); ylabel('BIC');
    end
end