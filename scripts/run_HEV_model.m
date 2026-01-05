function run_HEV_model()
    tic
    rng(123);

    % get folder containing this script
    this_file_path = mfilename('fullpath');
    this_folder = fileparts(this_file_path);
    repo_root = fullfile(this_folder, '..'); % move up folder up to repo main
    
    % create 'data' folder inside repo root if it doesn't exist
    data_dir = fullfile(repo_root, 'data');
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end

    % load or generate the LHS sample for parameters
    param_file = fullfile(data_dir, 'param_samples.mat');
    num_runs = 500;
    num_params = 8;  
    
    if exist(param_file, 'file')
        load(param_file, 'param_samples'); % param_samples is a 500x8 matrix
    else
        lhs_samples = lhsdesign(num_runs, num_params);

        % transform using inverse CDFs of desired distributions
        q_samples = 0.2 + (0.8 - 0.2)*lhs_samples(:,1); % [0.2, 0.8]
        phi_samples = 1 + (7 - 1)*lhs_samples(:,2); % all testing is done within 1 - 7 days
        epsilon_samples = 0.68 + (0.86 - 0.68) * lhs_samples(:,3); % [67.7 - 86.6]
        p_v_samples = 0.91 + (0.99 - 0.91) * lhs_samples(:,4); % [0.91 - 0.99]
        nu_max_samples  = 0.40 + (0.60 - 0.40) * lhs_samples(:,5); % [0.40, 0.60]
        wepsilon_max_samples = 0.40 + (0.80 - 0.40) * lhs_samples(:,6); % [0.40, 0.80]
        a_samples = 0.20 + (0.50 - 0.20) * lhs_samples(:,7); % [0.20, 0.50]
        c_samples = 1 + (5 - 1)*lhs_samples(:,8); % [1,5]
        
        % combine into matrix where each row is a run, columns are parameter vectors
        param_samples = [q_samples, phi_samples, epsilon_samples, p_v_samples, nu_max_samples, ...
                            wepsilon_max_samples, a_samples, c_samples];
        save(param_file, 'param_samples', '-v7.3');
    end

    % intervention selection
    disp('Interventions:');
    disp('1. Baseline');
    disp('2. Testing and Isolation')
    disp('3. Vaccination');
    disp('4. WASH (effect on theta)');
    disp('5. WASH (effect on delta)');
    disp('6. WASH (effect on both theta and delta)');
    disp('7. WASH (theta & delta) + Vaccination');
    intervention_choice = input('Enter the intervention number: ');
    if intervention_choice < 1 || intervention_choice > 7
        error('Intervention must be between 1 and 7.');
    end

    % parameters (defaults)
    theta = 0.0005;
    alpha = 0.9; % relative trans. of asymp. to symp. infection in water
    zeta = 0.9; % relative trans. of mild symp. to sev symp. infection in water
    delta = 0.1; % default clearance rate (will be replaced by sampling)
    sigma = 7.0/34.0;
    gamma_a = 7.0/38.0;
    gamma_m = 7.0/38.0;
    gamma_s = 7.0/38.0;
    c_default = 1.0;
    beta_h = 0.00547; % human-to-human transmissibility % when alpha=0.2, beta_h=0.045 and beta_w=0.099
    eta = 0.8; % relative trans. of asymp. to symp. (human-to-human)
    kappa = 0.9; % relative transmission of mild to severe symptomatic (human-to-human)
    beta_w = 0.026; % water-to-human transmissibility
    p = 0.8; % proportion asymptomatic
    p_v_default = 0.96; % proportion asymptomatic post vaccine
    % q = 0.7; % proportion isolated
    % phi = 0.5; % testing rate
    epsilon_sample = 1.0; % vaccine efficacy (default no vaccination)
    D = 0.0; % vaccination rate (default no vaccination)
    m = 0.85; % proportion mild symptomatic
    rho = 4.0; % time since vaccination (12 weeks)
    tau = 2.0/7.0; % testing delay (2 days in weeks)
    a_default = 0.2; % logistic growth rate for WASH_theta
    wepsilon_max_default = 0.7; % max WASH_theta effectiveness
    N = 10000; % total population
    v_max_default = 0.0; % max 60% can be vaccinated

    % map intervention choice to parameters
    [D, applyWASHontheta, applyWASHondelta] = pick_intervention(intervention_choice, D);

    % ICs
    W0 = 0.0; E0 = 10.0; S0 = N - E0; E_sv0 = 0.0; A0 = 0.0; A_sv0 = 0.0;
    I_m0 = 0.0; I_mv0 = 0.0; I_s0 = 0.0; I_sv0 = 0.0; Q0 = 0.0; V0 = 0.0; 
    E_v0 = 0.0; R0 = 0.0; R1_0 = 0.0; R_pv0 = 0.0; SInc0 = 0.0; VInc0 = 0.0;

    % reproduction number using default parameters
    R_0 = ( beta_h * (((p*eta)/gamma_a) + (m*(1-p)*kappa/gamma_m) + ((m-1)*(p-1)/gamma_s)) ) + ...
          ( (S0*beta_w*theta)/delta ) * ( (p*alpha)/gamma_a + (m*(1-p)*zeta/gamma_m) + ((m-1)*(p-1)/gamma_s) );
    R_h = beta_h * (((p*eta)/gamma_a) + (m*(1-p)*kappa/gamma_m) + ((m-1)*(p-1)/gamma_s));
    R_e = ((S0*beta_w*theta)/delta) * ((p*alpha)/gamma_a + (m*(1-p)*zeta/gamma_m) + ((m-1)*(p-1)/gamma_s));
    fprintf('Estimated R_0 (default) = %.3f\n', R_0);
    fprintf('Estimated Contribution to R_0 from Human Transmission (R_h) = %.3f\n', R_h);
    fprintf('Estimated Contribution to R_0 from Environmental Transmission (R_e) = %.3f\n', R_e);

    % set up time (weeks)
    dt = 1.0/8.0;
    t_max = 208.0;
    num_steps = round(t_max/dt);
    time = linspace(0, t_max, num_steps);

    % thresholds and simulation replications
    threshold_values = [1, 10, 25];
    num_thresholds = length(threshold_values);
    num_interventions = 7; 

    % aggregate simulation outputs (full dynamics and inc data) for current intervention
    % create cell array of replicates for each threshold
    local_results = cell(num_thresholds, 1);
    for t_idx = 1:num_thresholds
        local_results{t_idx} = cell(num_runs, 1);
    end

    % simulation loop (for thresholds and current intervention)
    for run_idx = 1:num_runs
        % sample sensitivity parameters
        q_sample = param_samples(run_idx, 1);
        phi_sample = param_samples(run_idx, 2);
        epsilon_sample = param_samples(run_idx, 3);
        p_v_sample = param_samples(run_idx, 4);
        nu_max_sample  = param_samples(run_idx, 5);
        wepsilon_max_sample = param_samples(run_idx, 6);
        a_sample = param_samples(run_idx, 7);
        c_sample = param_samples(run_idx, 8);

        % define which parameters are used according to intervention_choice:
        switch intervention_choice
            case 1 % baseline (no interventions)
                q_use = 0;
                phi_use = 0;
                epsilon_use = 1.0;
                p_v_use = p_v_default;
                v_max_use = 0; % no vaccination
                a_use = 0;
                wepsilon_max_use = 0;
                c_use = c_default;
            case 2 % testing/isolation: sample epsilon and nu_max
                q_use = q_sample;
                phi_use = phi_sample;
                epsilon_use = 1.0;
                p_v_use = p_v_default;
                v_max_use = 0;
                a_use = 0;
                wepsilon_max_use = 0;
                c_use = c_default;
            case 3 % vaccination only: sample epsilon, p_v, and nu_max
                q_use = q_sample;
                phi_use = phi_sample;
                epsilon_use = epsilon_sample;
                p_v_use = p_v_sample;
                v_max_use = nu_max_sample * N;
                a_use = 0;
                wepsilon_max_use = 0;
                c_use = c_default;
            case 4 % WASH-theta: sample a and wepsilon_max
                q_use = q_sample;
                phi_use = phi_sample;
                epsilon_use = 1.0;
                p_v_use = p_v_default; 
                v_max_use = 0;
                a_use = a_sample;
                wepsilon_max_use = wepsilon_max_sample;
                c_use = c_default;
            case 5 % WASH-delta: sample c
                q_use = q_sample;
                phi_use = phi_sample;
                epsilon_use = 1.0;
                p_v_use = p_v_default;
                v_max_use = 0;
                a_use = 0;
                wepsilon_max_use = 0;
                c_use = c_sample;
            case 6 % WASH-both: sample a, wepsilon_max, and c
                q_use = q_sample;
                phi_use = phi_sample;
                epsilon_use = 1.0;
                p_v_use = p_v_default;
                v_max_use = 0;
                a_use = a_sample;
                wepsilon_max_use = wepsilon_max_sample;
                c_use = c_sample;
            case 7 % WASH+Vacc: sample epsilon, nu_max, a, wepsilon_max, and c
                q_use = q_sample;
                phi_use = phi_sample;
                epsilon_use = epsilon_sample;
                p_v_use = p_v_sample;
                v_max_use = nu_max_sample * N;
                a_use = a_sample;
                wepsilon_max_use = wepsilon_max_sample;
                c_use = c_sample;
            otherwise
                error('Invalid intervention choice');
        end

        for t_idx = 1:num_thresholds
            T = threshold_values(t_idx);
            trigger_time_run = NaN;

            % preallocate arrays
            W_vec = zeros(1,num_steps); S_vec = zeros(1,num_steps); E_vec = zeros(1,num_steps); 
            E_sv_vec = zeros(1,num_steps); A_vec = zeros(1,num_steps); A_sv_vec = zeros(1,num_steps);
            I_m_vec = zeros(1,num_steps); I_mv_vec = zeros(1, num_steps); I_s_vec = zeros(1,num_steps); 
            I_sv_vec = zeros(1, num_steps); Q_vec = zeros(1,num_steps); V_vec = zeros(1,num_steps); 
            E_v_vec = zeros(1,num_steps); R_vec = zeros(1,num_steps); R1_vec = zeros(1,num_steps); 
            Rpv_vec = zeros(1,num_steps); nu_vec = zeros(1,num_steps); nu_E_vec = zeros(1,num_steps);
            nu_A_vec = zeros(1,num_steps); nu_R1_vec = zeros(1,num_steps); SInc_vec = zeros(1,num_steps); 
            VInc_vec = zeros(1,num_steps); vacc_vec_total = zeros(1,num_steps); vacc_vec_wasted = zeros(1, num_steps);
            wepsilon_vec = zeros(1,num_steps); wasted_vac_vec = zeros(1,num_steps);
            weekly_Inc = zeros(1, t_max); weekly_IncIs = zeros(1, t_max); weekly_IncIsQ = zeros(1, t_max);

            % apply ICs
            W_vec(1) = W0; E_vec(1) = E0; S_vec(1) = S0;
            E_sv_vec(1) = E_sv0; A_vec(1) = A0; A_sv_vec(1) = A_sv0;
            I_m_vec(1) = I_m0; I_mv_vec(1) = I_mv0; I_s_vec(1) = I_s0; 
            I_sv_vec(1) = I_sv0; Q_vec(1) = Q0; V_vec(1) = V0; 
            E_v_vec(1) = E_v0; R_vec(1) = R0; R1_vec(1) = R1_0; 
            Rpv_vec(1) = R_pv0; SInc_vec(1) = SInc0; VInc_vec(1) = VInc0; 
            vacc_vec_total(1) = 0; vacc_vec_wasted(1) = 0;
            t0 = []; % intialize WASH inflection time

            % main time loop
            for k = 2:num_steps
                W_n = W_vec(k-1); S_n = S_vec(k-1); E_n = E_vec(k-1);
                E_sv_n = E_sv_vec(k-1); A_n = A_vec(k-1); A_sv_n = A_sv_vec(k-1); 
                I_m_n = I_m_vec(k-1); I_mv_n = I_mv_vec(k-1); I_s_n = I_s_vec(k-1); 
                I_sv_n = I_sv_vec(k-1); Q_n = Q_vec(k-1); V_n = V_vec(k-1); 
                E_v_n = E_v_vec(k-1); R_n = R_vec(k-1); R1_n = R1_vec(k-1); 
                Rpv_n = Rpv_vec(k-1); vacc_total_n = vacc_vec_total(k-1); vacc_wasted_n = vacc_vec_wasted(k-1);

                delayed_index = max(1, (1:k) - round(tau/dt));
                I_s_delayed = I_s_vec(delayed_index);
                I_sv_delayed = I_sv_vec(delayed_index);
                integral_condition = q_use * phi_use * sum(I_s_delayed + I_sv_delayed) * dt; 

                % use real cumulative count (no delay) to log when T was actually crossed
                I_s_undelayed = I_s_vec(1:k);
                I_sv_undelayed = I_sv_vec(1:k);
                true_cumulative = q_use * phi_use * sum(I_s_undelayed + I_sv_undelayed) * dt;
                if isnan(trigger_time_run) && true_cumulative >= T
                    trigger_time_run = time(k);
                end

                if applyWASHontheta
                    if isempty(t0) && integral_condition >= T
                        t0 = time(k) + 2;
                    end
                
                    if ~isempty(t0) && time(k) >= t0
                        wepsilon_vec(k) = wepsilon_max_use / (1 + exp(-a_use*(time(k)-t0)));
                    else
                        wepsilon_vec(k) = 0;
                    end
                
                    % skip vaccination during WASH only cases:
                    if D == 0.0
                        [nu, nu_E, nu_A, nu_R1] = deal(0,0,0,0);
                    end
                else
                    wepsilon_vec(k) = 0;
                end

                % effect on delta: when true we scale delta
                if applyWASHondelta
                    if integral_condition >= T
                        current_delta = c_use * delta;
                    else
                        current_delta = delta;
                    end
                else
                    current_delta = delta;
                end

                % vaccination (when D neq 0), else zero
                [nu, nu_E, nu_A, nu_R1] = vaccination_flow(D, T, k, dt, tau, ...
                    I_s_vec, I_sv_vec, q_use, phi_use, gamma_s, S_n, E_n, A_n, R1_n, vacc_total_n, v_max_use);
    
                % implicit Euler updates
                theta_eff = theta * (1 - wepsilon_vec(k));
                shedding = theta_eff * (alpha*A_n + alpha*A_sv_n + zeta*I_m_n + zeta*I_mv_n + I_s_n + I_sv_n);
                force = beta_h * ((eta*A_n + eta*A_sv_n + kappa*I_m_n + kappa*I_mv_n + I_s_n + I_sv_n) / N) + beta_w * W_n;
    
                % compute u(t, rho)
                delayed_index_rho = max(1, k - round(rho/dt));
                current_integral = 0;
                for idx_R = delayed_index_rho:(k-1)
                    force_idx = beta_h * ((eta*A_vec(idx_R) + eta*A_sv_vec(idx_R) + kappa*I_m_vec(idx_R) + kappa*I_mv_vec(idx_R) + I_s_vec(idx_R) + I_sv_vec(idx_R)) / N) + beta_w * W_vec(idx_R);
                    current_integral = current_integral + force_idx * dt;
                end
                u_t_rho = nu_vec(delayed_index_rho) * exp(-current_integral);
    
                % difference equations
                W_np1 = (W_n + dt*shedding) / ( 1 + dt*(shedding + current_delta) ); 
                S_np1 = (S_n - dt*u_t_rho) / (1 + dt*force);
                E_np1 = (E_n + dt*(force*S_n - nu_E)) / (1 + dt*sigma);
                E_sv_np1 = (E_sv_n + dt*nu_E) / (1 + dt*sigma);
                A_np1 = (A_n + dt*(p*sigma*E_n - nu_A)) / (1 + dt*gamma_a);
                A_sv_np1 = (A_sv_n + dt*(nu_A + p*sigma*E_sv_n + p_v_use*sigma*E_v_n)) / (1 + dt*gamma_a);
                I_m_np1 = (I_m_n + dt*m*(1-p)*sigma*E_n) / (1 + dt*gamma_m);
                I_mv_np1 = (I_mv_n + dt*m*sigma* ((1-p)*E_sv_n + (1-p_v_use)*E_v_n) ) / (1+ dt*gamma_m);
                delayed_index_tau = max(1, k - round(tau/dt));
                I_s_delayed = I_s_vec(delayed_index_tau)*exp(-gamma_s*tau);
                I_sv_delayed = I_sv_vec(delayed_index_tau)*exp(-gamma_s*tau);
                I_s_np1 = (I_s_n + dt*(((1-m)*(1-p)*sigma*E_n) - (q_use*phi_use*I_s_delayed))) / (1 + dt*gamma_s);
                I_sv_np1 = (I_sv_n + dt*(((1-m)*(1-p)*sigma*E_sv_n) + ((1-m)*(1-p_v_use)*sigma*E_v_n) - (q_use*I_sv_delayed))) / (1 + dt*gamma_s);
                Q_np1 = (Q_n + dt*q_use*phi_use*exp(-gamma_s*tau)*(I_s_vec(delayed_index_tau) + I_sv_vec(delayed_index_tau))) / (1 + dt*gamma_s);
                V_np1 = (V_n + dt*u_t_rho) / (1 + dt*(1-epsilon_use)*force);
                E_v_np1 = (E_v_n + dt*(1-epsilon_use)*force*V_n) / (1 + dt*sigma);
                R_np1 = R_n + dt*(gamma_a*A_sv_n + gamma_m*I_mv_n + gamma_s*I_sv_n + gamma_s*Q_n);
                R1_np1 = R1_n + dt*(gamma_a*A_n + gamma_m*I_m_n + gamma_s*I_s_n - nu_R1);
                Rpv_np1 = Rpv_n + dt*nu_R1;
    
                SInc_np1 = dt*(force * S_n);
                VInc_np1 = dt*(1-epsilon_use)*force * V_n;
                vacc_total_np1 = vacc_total_n + dt*(nu + nu_E + nu_A + nu_R1);
                vacc_wasted_np1 = vacc_wasted_n + dt*(nu_E + nu_A + nu_R1);
    
                % store
                W_vec(k) = W_np1; S_vec(k) = S_np1; E_vec(k) = E_np1;
                E_sv_vec(k) = E_sv_np1; A_vec(k) = A_np1; A_sv_vec(k) = A_sv_np1;
                I_m_vec(k) = I_m_np1; I_mv_vec(k) = I_mv_np1; I_s_vec(k) = I_s_np1; 
                I_sv_vec(k) = I_sv_np1; Q_vec(k) = Q_np1; V_vec(k) = V_np1; 
                E_v_vec(k) = E_v_np1; R_vec(k) = R_np1; R1_vec(k) = R1_np1; 
                Rpv_vec(k) = Rpv_np1; nu_vec(k) = nu; nu_E_vec(k) = nu_E; 
                nu_A_vec(k) = nu_A; nu_R1_vec(k) = nu_R1;
                SInc_vec(k) = SInc_np1; VInc_vec(k) = VInc_np1; vacc_vec_total(k) = vacc_total_np1;
                vacc_vec_wasted(k) = vacc_wasted_np1;

                % weekly incidence accumulation
                week_idx = min(floor(time(k)), t_max);
                if week_idx >= 1 && week_idx <= t_max
                    weekly_Inc(week_idx) = weekly_Inc(week_idx) + (SInc_np1 + VInc_np1);
                    weekly_IncIs(week_idx) = weekly_IncIs(week_idx) + dt*(1-m)*((1-p)*sigma*E_n + (1-p)*sigma*E_sv_n + (1-p_v_use)*sigma*E_v_n);
                    weekly_IncIsQ(week_idx) = weekly_IncIsQ(week_idx) + dt*q_use*phi_use*exp(-gamma_s*tau)*(I_s_vec(delayed_index_tau) + I_sv_vec(delayed_index_tau));
                end
            end % end time loop
            
            if vacc_vec_total(end) > 0
                wasted_fraction = vacc_vec_wasted(end) / vacc_vec_total(end);
            else
                wasted_fraction = 0;
            end

            % cumulative sums
            Cum_Inc_vec = cumsum(weekly_Inc);
            Cum_IncIs_vec = cumsum(weekly_IncIs);
            Cum_IncIsQ_vec = cumsum(weekly_IncIsQ);
    
            % store evertthing in struct
            sim_result = struct('T', T, 'time', time, 'trigger_time', trigger_time_run, 'W', W_vec, ...
                'S', S_vec, 'E', E_vec, 'A', A_vec, 'A_sv', A_sv_vec, 'I_m', I_m_vec, 'I_mv', I_mv_vec,...
                'I_s', I_s_vec, 'I_sv', I_sv_vec, 'Q', Q_vec, 'V', V_vec, 'E_v', E_v_vec, ...
                'R', R_vec, 'R1', R1_vec, 'Rpv', Rpv_vec, 'vacc_total', vacc_vec_total, 'vacc_wasted', vacc_vec_wasted, ...
                'wasted_fraction', wasted_fraction, 'weekly_Inc', weekly_Inc, 'weekly_IncIs', weekly_IncIs, 'weekly_IncIsQ', weekly_IncIsQ, ...
                'Cum_Inc', Cum_Inc_vec, 'Cum_IncIs', Cum_IncIs_vec, 'Cum_IncIsQ', Cum_IncIsQ_vec);
    
            local_results{t_idx}{run_idx} = sim_result;
        end
    end % end simulation loop

    % aggregate average results over runs for each threshold
    weekly_time = 1:t_max;
    agg_results = cell(num_thresholds, 1);
    for t_idx = 1:num_thresholds
        % preallocate arrays, num steps is full time series
        avg_W = zeros(1, num_steps); avg_S = zeros(1, num_steps); avg_E = zeros(1, num_steps);
        avg_A = zeros(1, num_steps); avg_A_sv = zeros(1, num_steps); avg_I_m = zeros(1, num_steps);
        avg_I_mv = zeros(1, num_steps); avg_I_s = zeros(1, num_steps); avg_I_sv = zeros(1, num_steps);
        avg_Q = zeros(1, num_steps); avg_V = zeros(1, num_steps); avg_R = zeros(1, num_steps);
        avg_R1 = zeros(1, num_steps); avg_Rpv = zeros(1, num_steps);
        % for weekly and cumulative incidences (length = t_max)
        avg_weekly_Inc = zeros(1, t_max); avg_weekly_IncIs = zeros(1, t_max); avg_weekly_IncIsQ = zeros(1, t_max);
        avg_Cum_Inc = zeros(1, t_max); avg_Cum_IncIs = zeros(1, t_max); avg_Cum_IncIsQ = zeros(1, t_max);
        all_triggers = nan(num_runs,1); avg_vacc_total = zeros(1, num_steps); avg_vacc_wasted = zeros(1, num_steps);
        avg_wasted_fraction = 0.0;

        % loop over runs:
        for run_idx = 1:num_runs
            res = local_results{t_idx}{run_idx};
            avg_W = avg_W + res.W; avg_S = avg_S + res.S; avg_E = avg_E + res.E;
            avg_A = avg_A + res.A; avg_A_sv = avg_A_sv + res.A_sv;
            avg_I_m = avg_I_m + res.I_m; avg_I_mv = avg_I_mv + res.I_mv;
            avg_I_s = avg_I_s + res.I_s; avg_I_sv = avg_I_sv + res.I_sv;
            avg_Q = avg_Q + res.Q; avg_V = avg_V + res.V; avg_R = avg_R + res.R;
            avg_R1 = avg_R1 + res.R1; avg_Rpv = avg_Rpv + res.Rpv;
            avg_weekly_Inc = avg_weekly_Inc + res.weekly_Inc;
            avg_weekly_IncIs = avg_weekly_IncIs + res.weekly_IncIs;
            avg_weekly_IncIsQ = avg_weekly_IncIsQ + res.weekly_IncIsQ;
            avg_Cum_Inc = avg_Cum_Inc + res.Cum_Inc;
            avg_Cum_IncIs = avg_Cum_IncIs + res.Cum_IncIs;
            avg_Cum_IncIsQ = avg_Cum_IncIsQ + res.Cum_IncIsQ;
            all_triggers(run_idx) = res.trigger_time;
            avg_vacc_total = avg_vacc_total + res.vacc_total;
            avg_vacc_wasted = avg_vacc_wasted + res.vacc_wasted;
            avg_wasted_fraction = avg_wasted_fraction + res.wasted_fraction;
        end

        % divide by number of runs:
        avg_W = avg_W / num_runs; avg_S = avg_S / num_runs;
        avg_E = avg_E / num_runs; avg_A = avg_A / num_runs;
        avg_A_sv = avg_A_sv / num_runs; avg_I_m = avg_I_m / num_runs;
        avg_I_mv = avg_I_mv / num_runs; avg_I_s = avg_I_s / num_runs;
        avg_I_sv = avg_I_sv / num_runs; avg_Q = avg_Q / num_runs;
        avg_V = avg_V / num_runs; avg_R = avg_R / num_runs;
        avg_R1 = avg_R1 / num_runs; avg_Rpv = avg_Rpv / num_runs;
        avg_weekly_Inc = avg_weekly_Inc / num_runs;
        avg_weekly_IncIs = avg_weekly_IncIs / num_runs;
        avg_weekly_IncIsQ = avg_weekly_IncIsQ / num_runs;
        avg_Cum_Inc = avg_Cum_Inc / num_runs;
        avg_Cum_IncIs = avg_Cum_IncIs / num_runs;
        avg_Cum_IncIsQ = avg_Cum_IncIsQ / num_runs;
        avg_trigger_time = mean(all_triggers, 'omitnan');
        avg_vacc_total = avg_vacc_total / num_runs;
        avg_vacc_wasted = avg_vacc_wasted / num_runs;
        avg_wasted_fraction = avg_wasted_fraction / num_runs;
        
        % compute final summary values for interventions involving vaccination
        if ismember(intervention_choice, [3, 7])
            total_vacc_given = avg_vacc_total(end);
            total_vacc_wasted = avg_vacc_wasted(end);
            wasted_frac = total_vacc_wasted / total_vacc_given;
        else
            total_vacc_given = NaN;
            total_vacc_wasted = NaN;
            wasted_frac = NaN;
        end

        agg_results{t_idx} = struct('T', threshold_values(t_idx), 'time_full', time, 'W', avg_W, 'S', avg_S, 'E', avg_E, 'A', avg_A, ...
            'A_sv', avg_A_sv, 'I_m', avg_I_m, 'I_mv', avg_I_mv, 'I_s', avg_I_s, 'I_sv', avg_I_sv, 'Q', avg_Q, 'V', avg_V, 'R', avg_R, ...
            'R1', avg_R1, 'Rpv', avg_Rpv, 'weekly_Inc', avg_weekly_Inc, 'weekly_IncIs', avg_weekly_IncIs, 'weekly_IncIsQ', avg_weekly_IncIsQ, ...
            'Cum_Inc', avg_Cum_Inc, 'Cum_IncIs', avg_Cum_IncIs, 'Cum_IncIsQ', avg_Cum_IncIsQ, 'avg_trigger_time', avg_trigger_time, ...
            'total_vacc_given', total_vacc_given, 'total_vacc_wasted', total_vacc_wasted, 'wasted_fraction', wasted_frac);
    end

    % save to global .mat file (3 thresholds x 7 interventions)
    this_file_path = mfilename('fullpath');
    this_folder = fileparts(this_file_path);
    repo_root = fullfile(this_folder, '..');

    data_dir = fullfile(repo_root, 'data');
    if ~exist(data_dir, 'dir')
        mkdir(data_dir);
    end
    
    save_filename = fullfile(data_dir, 'HEV_results_all.mat');

    if exist(save_filename, 'file')
        load(save_filename, 'all_results', 'all_local_results', 'all_thresholds', 'all_interventions');
        if ~isequal(size(all_results), [num_thresholds, num_interventions])
            error('all_results dimension mismatch. Should be %dx%d.', num_thresholds, num_interventions);
        end
    else
        all_results = cell(num_thresholds, num_interventions);
        all_local_results = cell(num_thresholds, num_interventions);
        all_thresholds = threshold_values;
        all_interventions = {'Baseline', 'Testing/Isolation','Vaccination','WASH-theta','WASH-delta','WASH-both','WASH+Vacc'};
    end
    
    for t_idx = 1:num_thresholds
        all_results{t_idx, intervention_choice} = agg_results{t_idx};
        all_local_results{t_idx, intervention_choice} = local_results{t_idx}; % save per-run data
    end
    
    save(save_filename, 'all_results', 'all_local_results', 'all_thresholds', 'all_interventions', '-v7.3');
    fprintf('Aggregated and per-run data saved to: %s\n', save_filename);

    %%%%%%% PLOTTING DYNAMICS AND INCIDENCE
    % 1. model dynamics for an exmaple threshold
    selected_T = 1;
    sel_idx = find(cellfun(@(x) x.T == selected_T, agg_results), 1);
    if isempty(sel_idx)
        error('Threshold %d not found in aggregated results.', selected_T);
    end
    dyn = agg_results{sel_idx}; % dynamics for selected threshold

    figure;
    plot(dyn.time_full, dyn.S, 'b', 'LineWidth', 2, 'DisplayName','Susceptible'); hold on;
    plot(dyn.time_full, dyn.E, 'r', 'LineWidth', 2, 'DisplayName','Exposed');
    plot(dyn.time_full, dyn.A, 'g', 'LineWidth', 2, 'DisplayName','Asymptomatic');
    plot(dyn.time_full, dyn.A_sv, 'Color',[0.85 0.33 0.10], 'LineWidth', 2, 'DisplayName','Exposed (V)');
    plot(dyn.time_full, dyn.I_m, 'm', 'LineWidth', 2, 'DisplayName','Mild Inf');
    plot(dyn.time_full, dyn.I_mv, 'y', 'LineWidth', 2, 'DisplayName', 'Mild/Vacc');
    plot(dyn.time_full, dyn.I_s, 'k', 'LineWidth', 2, 'DisplayName','Severe Inf');
    plot(dyn.time_full, dyn.I_mv, 'Color', [0.72 0.50 0.72], 'LineWidth', 2, 'DisplayName', 'Severe/Vacc');
    plot(dyn.time_full, dyn.Q, 'c', 'LineWidth', 2, 'DisplayName','Isolated');
    plot(dyn.time_full, dyn.V, 'Color',[0.20 0.78 0.72], 'LineWidth', 2, 'DisplayName','Vaccinated');
    plot(dyn.time_full, dyn.R, 'Color',[0.93 0.69 0.13], 'LineWidth', 2, 'DisplayName','Recovered (Vacc)');
    plot(dyn.time_full, dyn.R1, 'Color',[0.00 0.50 0.50], 'LineWidth', 2, 'DisplayName','Recovered (Not Vacc)');
    plot(dyn.time_full, dyn.Rpv, 'Color',[0.50 0.00 0.50], 'LineWidth', 2, 'DisplayName','Recovered (Vacc)');
    xlabel('Time (weeks)'); ylabel('Population Count');
    title(sprintf('Averaged Dynamics for T = %d, Intervention: %s', selected_T, all_interventions{intervention_choice}));
    legend('show','Location','best');
    hold off;

    % 2. weekly incidence (overall, severe, identified Severe)
    blue_shades = [0,0,0.5; 0,0,1; 0.5,0.3,0.85; 0.20,0.78,0.72; 0.53,0.81,0.92];
    
    % overall weekly incidence: one curve per threshold
    figure; hold on;
    for t_idx = 1:num_thresholds
        plot(1:t_max, all_results{t_idx, intervention_choice}.weekly_Inc, 'Color', blue_shades(t_idx,:), 'LineWidth', 2, ...
             'DisplayName', ['T = ', num2str(all_results{t_idx, intervention_choice}.T)]);
    end
    xlabel('Time (weeks)'); ylabel('Weekly Incidence (Overall)');
    title(['Weekly Incidence for Intervention: ', all_interventions{intervention_choice}]);
    legend('show');
    hold off;

    % severe weekly incidence
    figure; hold on;
    for t_idx = 1:num_thresholds
        plot(1:t_max, all_results{t_idx, intervention_choice}.weekly_IncIs, 'Color', blue_shades(t_idx,:), 'LineWidth', 2, ...
             'DisplayName', ['T = ', num2str(all_results{t_idx, intervention_choice}.T)]);
    end
    xlabel('Time (weeks)'); ylabel('Weekly Incidence (Severe)');
    title(['Weekly Severe Incidence for Intervention: ', all_interventions{intervention_choice}]);
    legend('show');
    hold off;

    % identified severe weekly incidence (will be 0 for baseline)
    figure; hold on;
    for t_idx = 1:num_thresholds
        plot(1:t_max, all_results{t_idx, intervention_choice}.weekly_IncIsQ, 'Color', blue_shades(t_idx,:), 'LineWidth', 2, ...
             'DisplayName', ['T = ', num2str(all_results{t_idx, intervention_choice}.T)]);
    end
    xlabel('Time (weeks)'); ylabel('Weekly Incidence (Identified Severe)');
    title(['Weekly Identified Severe Incidence for Intervention: ', all_interventions{intervention_choice}]);
    legend('show');
    hold off;

    % 3. cumulative incidence plots (overall, severe, identified severe (Q is 0 for baseline))
    figure; hold on;
    for t_idx = 1:num_thresholds
        plot(1:t_max, all_results{t_idx, intervention_choice}.Cum_Inc, 'Color', blue_shades(t_idx,:), 'LineWidth', 2, ...
             'DisplayName', ['T = ', num2str(all_results{t_idx, intervention_choice}.T)]);
    end
    xlabel('Time (weeks)'); ylabel('Cumulative Incidence (Overall)');
    title(['Cumulative Incidence for Intervention: ', all_interventions{intervention_choice}]);
    legend('show');
    hold off;

    figure; hold on;
    for t_idx = 1:num_thresholds
        plot(1:t_max, all_results{t_idx, intervention_choice}.Cum_IncIs, 'Color', blue_shades(t_idx,:), 'LineWidth', 2, ...
             'DisplayName', ['T = ', num2str(all_results{t_idx, intervention_choice}.T)]);
    end
    xlabel('Time (weeks)'); ylabel('Cumulative Incidence (Severe)');
    title(['Cumulative Severe Incidence for Intervention: ', all_interventions{intervention_choice}]);
    legend('show');
    hold off;

    figure; hold on;
    for t_idx = 1:num_thresholds
        plot(1:t_max, all_results{t_idx, intervention_choice}.Cum_IncIsQ, 'Color', blue_shades(t_idx,:), 'LineWidth', 2, ...
             'DisplayName', ['T = ', num2str(all_results{t_idx, intervention_choice}.T)]);
    end
    xlabel('Time (weeks)'); ylabel('Cumulative Incidence (Identified Severe)');
    title(['Cumulative Identified Severe Incidence for Intervention: ', all_interventions{intervention_choice}]);
    legend('show');
    hold off;

    % sensitivity data to .csv: create a table for current intervention
    q_samples = param_samples(:,1);
    phi_samples = param_samples(:,2);
    epsilon_samples = param_samples(:,3);
    p_v_samples = param_samples(:,4);
    nu_max_samples  = param_samples(:,5); 
    wepsilon_max_samples = param_samples(:,6);
    a_samples = param_samples(:,7);
    c_samples = param_samples(:,8); 

    % create the table of 8 parameter columns:
    T_params = table(q_samples, phi_samples, epsilon_samples, p_v_samples, nu_max_samples,...
                        wepsilon_max_samples, a_samples, c_samples, ...
                        'VariableNames', {'q', 'phi', 'epsilon', 'p_v', 'nu_max',... 
                        'wepsilon_max', 'a', 'c'});
    
    % write them only if the CSV file does not exist:
    csv_filename = fullfile(data_dir, 'sensitivity_HEV_data.csv');
    if ~isfile(csv_filename)
        writetable(T_params, csv_filename); % if csv doesn't exist, write the parameter columns
        fprintf('Created CSV file with parameter columns: %s\n', csv_filename);
    else
        fprintf('CSV file already exists. Not rewriting parameter columns.\n'); % csv does exist, skip rewriting parameter columns
    end
    toc
end

%%%%%%% HELPER FUNCTIONS
% Intervention logic
function [D, applyWASHontheta, applyWASHondelta] = pick_intervention(choice, D_in)
    D = D_in;
    applyWASHontheta = false;
    applyWASHondelta = false;
    switch choice
        case 1
            fprintf('Intervention 1: Baseline\n');
        case 2
            fprintf('Intervention 2: Testing and Isolation\n');
        case 3
            D = 500.0;
            fprintf('Intervention 3: Vaccination\n');
        case 4
            applyWASHontheta = true;
            fprintf('Intervention 4: (WASH) Sanitation\n');
        case 5
            applyWASHondelta = true;
            fprintf('Intervention 5: (WASH) Water Treatment\n');
        case 6
            applyWASHontheta = true;
            applyWASHondelta = true;
            fprintf('Intervention 6: WASH (Sanitation and Water Treatment)\n');
        case 7
            D = 500.0;
            applyWASHontheta = true;
            applyWASHondelta = true;
            fprintf('Intervention 7: WASH + Vaccination\n');
        otherwise
            error('Invalid intervention choice.');
    end
end

% Vaccination logic
function [nu, nu_E, nu_A, nu_R1] = vaccination_flow(D, T, k, dt, tau, I_s_vec, I_sv_vec, q, phi, gamma_s, S_n, E_n, A_n, R1_n, vacc_total_n, v_max)
    if D == 0.0
        nu = 0; nu_E = 0; nu_A = 0; nu_R1 = 0;
        return;
    end
    delayed_index_tau = max(1, k - round(tau/dt));
    I_s_delayed_array = I_s_vec(1:delayed_index_tau);
    I_sv_delayed_array = I_sv_vec(1:delayed_index_tau);
    integral_condition = q * phi * exp(-gamma_s*tau) * sum(I_s_delayed_array + I_sv_delayed_array) * dt;
    if vacc_total_n >= v_max
        nu = 0; nu_E = 0; nu_A = 0; nu_R1 = 0;
    else
        if integral_condition < T
            nu = 0; nu_E = 0; nu_A = 0; nu_R1 = 0;
        else
            if S_n + E_n + A_n + R1_n == 0
                nu = 0; nu_E = 0; nu_A = 0; nu_R1 = 0;
            else
                nu = min(S_n, (S_n*D)/(S_n + E_n + A_n + R1_n));
                nu_E = min(E_n, (E_n/S_n)*nu);
                nu_A = min(A_n, (A_n/S_n)*nu);
                nu_R1 = min(R1_n, (R1_n/S_n)*nu);
            end
        end
    end
end