function analyze_HEV_results()
    % set path and load data:
    this_file_path = mfilename('fullpath');
    this_folder = fileparts(this_file_path);
    repo_root = fullfile(this_folder, '..');

    data_dir = fullfile(repo_root, 'data');
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    mat_file = fullfile(data_dir, 'HEV_results_all.mat');
    load(mat_file, 'all_results', 'all_local_results', 'all_thresholds', 'all_interventions');
    num_thresholds = length(all_thresholds);

    % pick threshold
    disp('Available thresholds are:'); 
    disp(all_thresholds);
    chosen_T = input('Pick a threshold value (must be 1, 10, or 25): ');
    idx_T = find(all_thresholds == chosen_T, 1);
    if isempty(idx_T)
        error('Threshold %d not found in all_thresholds.', chosen_T);
    end

    % pick two interventions
    disp('Interventions:');
    for i = 1:numel(all_interventions)
        fprintf('  %d) %s\n', i, all_interventions{i});
    end
    baseline_col = input('Pick the "baseline" intervention column (1...7): ');
    if baseline_col < 1 || baseline_col > 7
        error('Column must be between 1 and 7.');
    end
    intervention_col = input('Pick the "intervention" column (1...7): ');
    if intervention_col < 1 || intervention_col > 7
        error('Column must be between 1 and 7.');
    end

    % retrieve aggregated results for the selected threshold and interventions
    baseline_res = all_results{idx_T, baseline_col}; 
    intervention_res = all_results{idx_T, intervention_col};
    
    Tval = baseline_res.T;  
    tvec = baseline_res.time_full; % common time vector
    
    % mean trigger time (in weeks) for this threshold
    t_base = baseline_res.avg_trigger_time;
    t_int  = intervention_res.avg_trigger_time;
    fprintf('Mean trigger time (baseline) = %.2f weeks\n', t_base);
    fprintf('Mean trigger time (intervention)= %.2f weeks\n\n', t_int);

    % find first time index where Cum_IncIsQ >= T
    ix_base = find(baseline_res.Cum_IncIsQ >= Tval, 1, 'first');
    ix_int  = find(intervention_res.Cum_IncIsQ >= Tval, 1, 'first');
    
    if isempty(ix_base)
        fprintf('Baseline never crossed T=%d\n', Tval);
    else
        fprintf('Baseline rollout (threshold) at t = %.2f weeks\n', tvec(ix_base));
    end
    
    if isempty(ix_int)
        fprintf('Intervention never crossed T=%d\n', Tval);
    else
        fprintf('Intervention rollout at t = %.2f weeks\n', tvec(ix_int));
    end

    % compute aggregated metrics (from the averaged curves)
    total_inf_baseline = baseline_res.Cum_Inc(end);
    total_inf_interv = intervention_res.Cum_Inc(end);
    tot_prevented = total_inf_baseline - total_inf_interv;
    total_sev_baseline = baseline_res.Cum_IncIs(end);
    total_sev_interv = intervention_res.Cum_IncIs(end);
    sev_prevented = total_sev_baseline - total_sev_interv;

    fprintf('(1) Total cases baseline: %.4f\n', total_inf_baseline);
    fprintf('(1) Total cases intervention: %.4f\n', total_inf_interv);
    fprintf('(1) Total cases prevented: %.4f\n', tot_prevented);
    fprintf('(1) Total severe cases baseline: %.4f\n', total_sev_baseline);
    fprintf('(1) Total severe cases intervention: %.4f\n', total_sev_interv);
    fprintf('(1) Severe cases prevented: %.4f\n', sev_prevented);
    
    rel_reduction_sev = (total_sev_baseline - total_sev_interv) / total_sev_baseline * 100;
    rel_reduction_tot = (total_inf_baseline - total_inf_interv) / total_inf_baseline * 100;
    fprintf('\n Results for T = %d \n', chosen_T);
    fprintf('Comparing: %s (col = %d) vs. %s (col = %d)\n', ...
        all_interventions{baseline_col}, baseline_col, all_interventions{intervention_col}, intervention_col);
    fprintf('\n(2) Relative Reduction in total cases: %.2f %%\n', rel_reduction_tot);
    fprintf('\n(2) Relative Reduction in severe cases: %.2f %%\n', rel_reduction_sev);
    
    % append aggregated plotting metrics to a .mat file
    metrics.T = chosen_T;
    metrics.baseline_intervention = baseline_col;
    metrics.intervention = intervention_col;
    metrics.trigger_base = t_base;
    metrics.trigger_interv = t_int;
    metrics.weekly_inc_severe_baseline = baseline_res.weekly_IncIs;
    metrics.weekly_inc_severe_intervention = intervention_res.weekly_IncIs;
    metrics.weekly_inc_overall_baseline = baseline_res.weekly_Inc;
    metrics.weekly_inc_overall_intervention = intervention_res.weekly_Inc;
    metrics.rel_reduction_sev = rel_reduction_sev;
    metrics.rel_reduction_tot = rel_reduction_tot;
    
    plot_file = fullfile(data_dir, 'HEV_results_plotting.mat');
    if exist(plot_file, 'file')
        load(plot_file, 'metrics_data');
    else
        metrics_data = [];
    end
    metrics_data = [metrics_data; metrics];
    save(plot_file, 'metrics_data', '-v7.3');
    fprintf('Aggregated metrics saved to: %s\n', plot_file);
    
    % compute the relative reduction (RR) per simulation run for the chosen threshold
    % all_local_results is a 3x7 cell array where each cell contains a cell array of run structs
    % use chosen threshold idx_T and baseline (col=1) for comparing with interventions (2-7)
    num_runs = length(all_local_results{idx_T,1}); % number of runs for the chosen threshold
    rr_vec = zeros(num_runs, 1);
    for i = 1:num_runs
        base_cum = all_local_results{idx_T, 1}{i}.Cum_Inc(end);
        interv_cum = all_local_results{idx_T, intervention_col}{i}.Cum_Inc(end);
        if base_cum > 0
            rr_vec(i) = (base_cum - interv_cum) / base_cum * 100;
        else
            rr_vec(i) = NaN;
        end
    end

    % create column for selected threshold and intervention pair
    RR_current = table(rr_vec, 'VariableNames', {sprintf('RR_T%d_I%d', chosen_T, intervention_col)});
    csv_filename = fullfile(data_dir, 'sensitivity_HEV_data.csv'); % append RR column to CSV file

    if isfile(csv_filename)
        sens_data = readtable(csv_filename);
        % remove any existing column for the current threshold and intervention to avoid duplicates
        col_to_remove = sprintf('RR_T%d_I%d', chosen_T, intervention_col);
        if ismember(col_to_remove, sens_data.Properties.VariableNames)
            sens_data.(col_to_remove) = [];
        end
        combined_data = [sens_data, RR_current];
    else
        combined_data = RR_current;
    end
    writetable(combined_data, csv_filename);
    fprintf('Sensitivity RR data appended to CSV: %s\n', csv_filename);

    % print vaccination summary only for vaccination scenarios (3 or 7)
    if ismember(intervention_col, [3, 7])
        total_vacc_given = intervention_res.total_vacc_given;
        total_vacc_wasted = intervention_res.total_vacc_wasted;
        wasted_fraction = intervention_res.wasted_fraction; % scalar
    
        fprintf('\n--- Vaccination Summary for %s (Intervention %d) ---\n', ...
            all_interventions{intervention_col}, intervention_col);
        fprintf('Total vaccines given: %.2f\n', total_vacc_given);
        fprintf('Total vaccines wasted: %.2f\n', total_vacc_wasted);
        fprintf('Fraction wasted: %.2f %%\n', wasted_fraction * 100);
    end
end