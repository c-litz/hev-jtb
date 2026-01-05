clear; clc;

% set path, load data, organize data:
this_file_path = mfilename('fullpath');
this_folder = fileparts(this_file_path);
repo_root = fullfile(this_folder, '..');

data_dir = fullfile(repo_root, 'data'); 
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

mat_file = fullfile(data_dir, 'HEV_results_plotting.mat');
load(mat_file, 'metrics_data');

all_T = unique([metrics_data.T]); % extract thresholds 
all_interv_temp = unique([metrics_data.intervention]); % [2 3 4 5 6 7]
if ~ismember(1, all_interv_temp)
    all_interv = [1, all_interv_temp]; % force in baseline as intervention 1
else
    all_interv = all_interv_temp;
end
num_T = numel(all_T);
num_int = numel(all_interv);

all_interv = unique([metrics_data.intervention]); % extract interventions 2-7

% preallocate matrices
RR_sev = nan(num_T, num_int-1); % RR severe incidence
RR_tot = nan(num_T, num_int-1); % RR total incidence
weekly_inc_severe = cell(num_int, num_T); % weekly severe incidence
weekly_inc_overall = cell(num_int, num_T); % weekly overall incidence
cum_inc_severe = cell(num_int, num_T); % cumulative severe inc
cum_inc_overall = cell(num_int, num_T); % cumulative overall inc

% loop through metrics_data and fill matrices
for i = 1:length(metrics_data)
    m_data = metrics_data(i);
    t_idx = find(all_T == m_data.T); % find corresponding indices for T and intervention
    int_idx = find(all_interv == m_data.intervention); % struct = numeric
    
    if ~isempty(t_idx) && ~isempty(int_idx)
        RR_sev(t_idx, int_idx) = m_data.rel_reduction_sev;
        RR_tot(t_idx, int_idx) = m_data.rel_reduction_tot;
    end
    % incidences for plotting
    if m_data.intervention == 1
        weekly_inc_severe{int_idx, t_idx} = m_data.weekly_inc_severe_baseline;
        weekly_inc_overall{int_idx, t_idx}  = m_data.weekly_inc_overall_baseline;
    else
        weekly_inc_severe{int_idx, t_idx} = m_data.weekly_inc_severe_intervention;
        weekly_inc_overall{int_idx, t_idx}  = m_data.weekly_inc_overall_intervention;
    end
    
    % cumulative incidence from weekly curves:
    cum_inc_severe{int_idx, t_idx} = cumsum(weekly_inc_severe{int_idx, t_idx});
    cum_inc_overall{int_idx, t_idx} = cumsum(weekly_inc_overall{int_idx, t_idx});
end

mat_file = fullfile(data_dir, 'HEV_results_all.mat');
load(mat_file, 'all_local_results', 'all_thresholds', 'all_interventions');
% all_local_results is a num_thresholds × num_interventions cell array, each cell is 1×num_runs cell array

num_T = numel(all_thresholds);
num_int = numel(all_interventions);

% compute trigger times:
avg_trigger_times = nan(num_T, num_int);
for j = 1:num_T
    Tval = all_thresholds(j);
    for m = 1:num_int
        runs = all_local_results{j, m}; % 1×num_runs cell array
        nRuns = numel(runs);
        times = nan(nRuns, 1);
        for r = 1:nRuns
            runStruct = runs{r};
            if isfield(runStruct, 'trigger_time') && ~isnan(runStruct.trigger_time)
                times(r) = runStruct.trigger_time;
            end
        end
        avg_trigger_times(j, m) = mean(times, 'omitnan');
    end
end

% output for trigger times: table with row=threshold & col=intervention name
T = array2table(avg_trigger_times, 'RowNames', compose("T=%d", all_thresholds), ...
    'VariableNames', all_interventions);
disp(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General plotting necessities:
intervs = {'Testing & Isolation','Vaccination','Sanitation','Water Treatment','WASH','Vaccination & WASH'};
blue_shades = [0,0,0.5; 0,0,1; 0.5,0.3,0.85; 0.20,0.78,0.72; 0.53,0.81,0.92; 0.5176,0.5176,0.8314];
ll = 1.001;
t_idx_T1 = find(all_thresholds == 1);
t_idx_T25 = find(all_thresholds == 25);
t = 1:208;
interv_colors = blue_shades;

% BASELINE DATA FOR PLOTTING:
load(mat_file, 'all_results', 'all_thresholds', 'all_interventions'); % same matfile as above
idx_T = find(all_thresholds == 25, 1); % baseline is indepdendent of T, any T is ok
baseline_col = 1;
baseline_res = all_results{idx_T, baseline_col};

weekly_inc_sev_baseline = baseline_res.weekly_IncIs;
cum_inc_sev_baseline = baseline_res.Cum_IncIs;
weekly_inc_overall_baseline = baseline_res.weekly_Inc;
cum_inc_overall_baseline = baseline_res.Cum_Inc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% GROUPED BAR CHARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tile = tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');
ylabel(tile, 'Relative Reduction (%)', 'FontSize', 18);

T_subset_idx = ismember(all_T, [1, 25]);
T_subset = all_T(T_subset_idx);
RR_sev_subset = RR_sev(T_subset_idx, :);
RR_tot_subset = RR_tot(T_subset_idx, :);

% total cases
nexttile(tile,1);
b1 = bar(T_subset, RR_tot_subset, 'grouped');
for k = 1:numel(b1)
    b1(k).FaceColor = blue_shades(k, :);
end

ylim([0 100]);
ax = gca;
ax.FontSize = 18; 
grid on;
ax.XGrid = 'off';
ax.YGrid = 'on';
set(gca, 'XTickLabel', []); % no x labels
box off;
text(0.01, 1.11, '(A)', 'Units', 'normalized', ...
     'FontWeight', 'bold', 'FontSize', 26, ...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

lgd = legend(b1, intervs, 'Location', 'northeast', 'Box', 'off', ...
    'FontSize', 12, 'NumColumns', 2);

% severe cases
nexttile(tile,2);
b2 = bar(T_subset, RR_sev_subset, 'grouped');
for k = 1:numel(b2)
    b2(k).FaceColor = blue_shades(k, :);
end

xlabel('Threshold T','FontSize', 14);
ylim([0 100]);
grid on;
ax = gca;
ax.FontSize = 18;
grid on;
ax.XGrid = 'off'; 
ax.YGrid = 'on';
box off;
text(0.01, 1.11, '(B)', 'Units', 'normalized', 'FontWeight', 'bold', ...
     'FontSize', 26, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

fig_dir = fullfile(repo_root, 'figures');  % or 'Figures' if you prefer capital F
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

outfile = fullfile(fig_dir, 'Figure4.pdf');
drawnow;
set(gcf, 'Units', 'inches', 'Position', [1 1 11 12]);
exportgraphics(gcf, outfile, 'ContentType', 'vector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% WATER CONTAMINATION INDEX AND SEVERE INCIDENCE %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_labels = ["Baseline", "Testing & Isolation", "Vaccination", ...
              "Sanitation", "Water Treatment", "WASH", "Vaccination & WASH"];
T_vals = [1, 25];
t_idx_vals = [t_idx_T1, t_idx_T25];

fig = figure('Units','inches','Position',[1 1 9.5 40]);
set(fig, 'Color', 'w');

% tile layout:
lefts = [0.11, 0.51];
width = 0.32;
top_margin = 0.07;
bottom_margin = 0.07;
v_gap = 0.025;
num_plots = 7;
total_gap = (num_plots - 1) * v_gap;
available_height = 1 - top_margin - bottom_margin - total_gap;
height = available_height / num_plots;

for T_plot = 1:2
    t_idx = t_idx_vals(T_plot);
    left = lefts(T_plot);

    for int_idx = 1:num_plots
        bottom = bottom_margin + (num_plots - int_idx) * (height + v_gap);
        ax = axes(fig, 'Position', [left, bottom, width, height]);
        hold(ax, 'on');

        % data
        if int_idx == 1
            inc_curve = cum_inc_sev_baseline;
            WSI_curve = all_results{t_idx, 1}.W;
        else
            inc_curve = all_results{t_idx, int_idx}.Cum_IncIs;
            WSI_curve = all_results{t_idx, int_idx}.W;
        end
        N = 10000;
        inc_curve = (inc_curve / N) * 100;

        % need the water index weekly
        num_weeks = floor(length(WSI_curve)/8);
        WSI_downsampled = mean(reshape(WSI_curve(1:num_weeks*8), 8, num_weeks), 1);
        WSI_downsampled = min(max(WSI_downsampled, 0), 1);

        % background fill
        yl = [0 3];
        num_bins = 8;
        bin_edges = linspace(0, 1, num_bins + 1);
        cmap = jet(num_bins);
        for i = 1:length(t)
            wval = WSI_downsampled(i);
            [~, bin_idx] = histc(wval, bin_edges);
            bin_idx = min(max(bin_idx, 1), num_bins);
            col = cmap(bin_idx, :);
            fill([t(i)-0.5 t(i)+0.5 t(i)+0.5 t(i)-0.5], [yl(1) yl(1) yl(2) yl(2)], ...
                col, 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Parent', ax);
        end

        % plot incidence curve
        plot(ax, t, inc_curve, 'k-', 'LineWidth', 1.5);
        xlim(ax, [0 208]);
        ylim(ax, yl);
        yticks(ax, 0:1.5:3);
        yticklabels(ax, {'0', '1.5', '3'});
        set(ax, 'FontSize', 14);
        grid on;
        title(ax, int_labels{int_idx}, 'FontWeight','bold', 'FontSize', 12);

        % draw trigger line
        if int_idx > 1
            trigger_time = avg_trigger_times(t_idx, int_idx);

            if T_vals(T_plot) == 1
                trigger_color = [1 1 1];
            elseif T_vals(T_plot) == 25
                trigger_color = [0 0 0];
            else
                trigger_color = [0 0 0];
            end

            xline(ax, trigger_time, '--', 'Color', trigger_color, ...
                'LineWidth', 1.4, 'Label', sprintf('T=%d', T_vals(T_plot)), ...
                'LabelOrientation','horizontal', 'LabelVerticalAlignment','bottom', ...
                'FontWeight','bold', 'FontSize', 11);
        end

        % x-axis labels
        if int_idx == num_plots
            xlabel(ax, 'Time (Weeks)', 'FontSize', 18);
        else
            set(ax, 'XTickLabel', []);
        end
    end

    % (A) and (B) titles at top left
    annotation(fig, 'textbox', 'String', sprintf('(%s)', char(double('A') + T_plot - 1)), ...
        'Units', 'normalized', 'Position', [left, 0.945, 0.05, 0.03], ...
        'FontSize', 20, 'FontWeight','bold', 'EdgeColor','none');
end

% y-axis 
y_label_ax = axes('Position', [0.04, 0.35, 0.03, 0.3], 'Visible', 'off');
text(0.5, 0.5, 'Proportion of Population Infected with Severe Disease (%)', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 18, ...
    'Rotation', 90, 'Parent', y_label_ax);

% colorbar on far right
cb_ax = axes(fig, 'Position', [0.82, 0.15, 0.01, 0.68], 'Visible', 'off');
colormap(cb_ax, jet);
cb = colorbar(cb_ax, 'Orientation', 'vertical', 'Ticks', 0:0.1:1, 'FontSize', 14);
cb.Label.String = 'Water Contamination Index';
cb.Label.Rotation = 90;
cb.Label.FontSize = 18;

% export
fig_dir = fullfile(repo_root, 'figures');  % figures directory
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
outname = fullfile(fig_dir, sprintf('Figure3.pdf', T_vals(T_plot)));
exportgraphics(gcf, outname, 'ContentType', 'vector');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WATER CONTAMINATION INDEX AND OVERALL INCIDENCE %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_labels = ["Baseline", "Testing & Isolation", "Vaccination", ...
              "Sanitation", "Water Treatment", "WASH", "Vaccination & WASH"];
T_vals = [1, 25];
t_idx_vals = [t_idx_T1, t_idx_T25];

fig = figure('Units','inches','Position',[1 1 9.5 40]);
set(fig, 'Color', 'w');

% tile layout:
lefts = [0.11, 0.51];
width = 0.32;
top_margin = 0.07;
bottom_margin = 0.07;
v_gap = 0.025;
num_plots = 7;
total_gap = (num_plots - 1) * v_gap;
available_height = 1 - top_margin - bottom_margin - total_gap;
height = available_height / num_plots;

for T_plot = 1:2
    t_idx = t_idx_vals(T_plot);
    left = lefts(T_plot);

    for int_idx = 1:num_plots
        bottom = bottom_margin + (num_plots - int_idx) * (height + v_gap);
        ax = axes(fig, 'Position', [left, bottom, width, height]);
        hold(ax, 'on');

        % data
        if int_idx == 1
            inc_curve = cum_inc_overall_baseline;
            WSI_curve = all_results{t_idx, 1}.W;
        else
            inc_curve = all_results{t_idx, int_idx}.Cum_Inc;
            WSI_curve = all_results{t_idx, int_idx}.W;
        end
        N = 10000;
        inc_curve = (inc_curve / N) *100;

        % need the water index weekly
        num_weeks = floor(length(WSI_curve)/8);
        WSI_downsampled = mean(reshape(WSI_curve(1:num_weeks*8), 8, num_weeks), 1);
        WSI_downsampled = min(max(WSI_downsampled, 0), 1);

        % background fill
        yl = [0 100];
        num_bins = 8;
        bin_edges = linspace(0, 1, num_bins + 1);
        cmap = jet(num_bins);
        for i = 1:length(t)
            wval = WSI_downsampled(i);
            [~, bin_idx] = histc(wval, bin_edges);
            bin_idx = min(max(bin_idx, 1), num_bins);
            col = cmap(bin_idx, :);
            fill([t(i)-0.5 t(i)+0.5 t(i)+0.5 t(i)-0.5], [yl(1) yl(1) yl(2) yl(2)], ...
                col, 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Parent', ax);
        end

        % plot incidence curve
        plot(ax, t, inc_curve, 'k-', 'LineWidth', 1.5);
        xlim(ax, [0 208]);
        ylim(ax, [0 100]); 
        yticks(ax, [0 50 100]);
        yticklabels(ax, {'0', '50', '100'}); 
        set(ax, 'FontSize', 14);
        grid on;
        title(ax, int_labels{int_idx}, 'FontWeight','bold', 'FontSize', 12);

        % draw trigger line
        if int_idx > 1
            trigger_time = avg_trigger_times(t_idx, int_idx);

            if T_vals(T_plot) == 1
                trigger_color = [1 1 1];
            elseif T_vals(T_plot) == 25
                trigger_color = [0 0 0];
            else
                trigger_color = [0 0 0];
            end

            xline(ax, trigger_time, '--', 'Color', trigger_color, ...
                'LineWidth', 1.4, 'Label', sprintf('T=%d', T_vals(T_plot)), ...
                'LabelOrientation','horizontal', 'LabelVerticalAlignment','bottom', ...
                'FontWeight','bold', 'FontSize', 11);
        end

        % x-axis labels
        if int_idx == num_plots
            xlabel(ax, 'Time (Weeks)', 'FontSize', 18);
        else
            set(ax, 'XTickLabel', []);
        end
    end

    % (A) and (B) titles at top left
    annotation(fig, 'textbox', 'String', sprintf('(%s)', char(double('A') + T_plot - 1)), ...
        'Units', 'normalized', 'Position', [left, 0.945, 0.05, 0.03], ...
        'FontSize', 20, 'FontWeight','bold', 'EdgeColor','none');
end

% y-axis 
y_label_ax = axes('Position', [0.04, 0.35, 0.03, 0.3], 'Visible', 'off');
text(0.5, 0.5, 'Proportion of Population Infected (%)', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 18, ...
    'Rotation', 90, 'Parent', y_label_ax);

% colorbar
cb_ax = axes(fig, 'Position', [0.82, 0.15, 0.01, 0.68], 'Visible', 'off');
colormap(cb_ax, jet);
cb = colorbar(cb_ax, 'Orientation', 'vertical', 'Ticks', 0:0.1:1, 'FontSize', 14);
cb.Label.String = 'Water Contamination Index';
cb.Label.Rotation = 90;
cb.Label.FontSize = 18;

% export
fig_dir = fullfile(repo_root, 'figures');  % figures directory
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
outname = fullfile(fig_dir, sprintf('Figure2.pdf', T_vals(T_plot)));
exportgraphics(gcf, outname, 'ContentType', 'vector');
