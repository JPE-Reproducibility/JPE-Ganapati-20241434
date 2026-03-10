function figures_gen(baseline_path, output_path, version, qq, varargin)
    warning('off', 'MATLAB:legend:IgnoringExtraEntries');
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'specifications', {'LogNormal', 'MEstimates'}, @iscell);
    addParameter(p, 'save_separate', true, @islogical);
    addParameter(p, 'save_combined', true, @islogical);
    addParameter(p, 'plot_gft', false, @islogical);
    addParameter(p, 'gmm_estimators', {struct('basis','Spline','knots',5)}, @iscell);
    addParameter(p, 'ylim_theta', [0 8], @(x) isnumeric(x) && length(x)==2);
    addParameter(p, 'xlim_gft', [0 8], @(x) isnumeric(x) && length(x)==2);
    addParameter(p, 'gft_tick_delta', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'gft_binwidth', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'file_suffix', '', @ischar);
    addParameter(p, 'legend_labels', {}, @iscell);
    addParameter(p, 'legend_position', 'southeast', @ischar);
    addParameter(p, 'font_size', 12, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'legend_font_size', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    parse(p, varargin{:});
    
    specs = p.Results.specifications;
    save_separate = p.Results.save_separate;
    save_combined = p.Results.save_combined;
    plot_gft = p.Results.plot_gft;
    gmm_estimators = p.Results.gmm_estimators;
    ylim_theta = p.Results.ylim_theta;
    xlim_gft = p.Results.xlim_gft;
    gft_tick_delta = p.Results.gft_tick_delta;
    gft_binwidth = p.Results.gft_binwidth;
    file_suffix = p.Results.file_suffix;
    legend_labels = p.Results.legend_labels;
    legend_position = p.Results.legend_position;
    font_size = p.Results.font_size;
    legend_font_size = p.Results.legend_font_size;
    
    % If legend_font_size not specified, default to 1.5x the font_size
    if isempty(legend_font_size)
        legend_font_size = font_size * 1.5;
    end
    
    % Helper function to get qq specification
    function q = getq(spec)
        if isstruct(qq) && isfield(qq, spec)
            q = string(qq.(spec));
        else
            q = string(qq);
        end
    end
    
    % Specification mapping
    spec_map = struct(...
        'LogNormal', struct('label', 'Log-Normal (Decreasing Elasticity)', 'color', [0.00 0.60 0.50]), ...
        'Estimates', struct('label', 'Increasing Elasticity', 'color', [0.00 0.20 0.55]), ...
        'MPareto', struct('label', 'Modified Pareto', 'color', [0.93 0.69 0.13]), ...
        'MEstimates', struct('label', 'Modified Estimates', 'color', [0.49 0.18 0.56]));
    
    % Load and process data for each specification
    data = struct();
    for i = 1:length(specs)
        spec = specs{i};
        data.(spec) = load_specification_data(baseline_path, version, spec, ...
                                               getq(spec), gmm_estimators{1});
        
        % Load additional GMM estimators if more than one specified
        if length(gmm_estimators) > 1
            data.(spec).gmm_extra = load_gmm_estimators(baseline_path, version, ...
                spec, gmm_estimators(2:end));
        end
    end
    
    % Plot elasticity estimates (GMM vs QQ, with additional GMM if specified)
    if ~isempty(specs)
        plot_elasticity_estimates(data, specs, spec_map, output_path, ...
            save_separate, save_combined, gmm_estimators, ylim_theta, ...
            file_suffix, legend_labels, legend_position, font_size, legend_font_size);
    end
    
    % Plot GFT histograms if requested
    if plot_gft && ~isempty(specs)
        plot_gft_histograms(data, specs, spec_map, output_path, save_separate, xlim_gft, gft_tick_delta, gft_binwidth, font_size, legend_font_size);
    end
end

function data = load_specification_data(baseline_path, version, spec, q, gmm_baseline)
    % Load GMM baseline data - use the first GMM estimator as baseline
    baseline_file = fullfile(baseline_path, ...
        sprintf("Simulated_%s_Estimates_%s_%d_%s.mat", spec, ...
                gmm_baseline.basis, gmm_baseline.knots, version));
    bl = load(baseline_file);
    
    % Load QQ data
    qq_file = fullfile(baseline_path, ...
        sprintf("Simulated_%s_Estimates_QQLogNormal_1_%s_%s.mat", spec, q, version));
    qq = load(qq_file);
    
    % Compute support
    lnn = linspace(prctile(bl.Rn_mat, 0.5, 'all'), ...
                   prctile(bl.Rn_mat, 99.5, 'all'), 500)';
    
    % Extract medians
    median_bl = median(bl.p.theta_n, 2);
    median_qq = median(qq.p.theta_n, 2);
    median_tr = median(bl.p.theta_true, 2);
    
    % Compute confidence bands
    low = 2.5; high = 97.5;
    bl_LB = prctile(bl.p.theta_n, low, 2);
    bl_UB = prctile(bl.p.theta_n, high, 2);
    qq_LB = prctile(qq.p.theta_n, low, 2);
    qq_UB = prctile(qq.p.theta_n, high, 2);
    
    % Extract GFT data if available
    if isfield(bl.p, 'w_hat') && isfield(bl.p, 'w_hat_true')
        % Remove complex values before computing MSE
        w_hat_bl = bl.p.w_hat;
        w_hat_true_bl = bl.p.w_hat_true;
        
        % Filter out rows (countries) with non-zero imaginary parts
        tol = 1e-12;
        valid_bl = all(abs(imag(w_hat_bl)) <= tol, 2) & all(abs(imag(w_hat_true_bl)) <= tol, 2);
        w_hat_bl = real(w_hat_bl(valid_bl, :));
        w_hat_true_bl = real(w_hat_true_bl(valid_bl, :));
        
        % Compute MSE (mean over countries for each simulation)
        gft_bl_mse = mean((w_hat_bl - w_hat_true_bl).^2, 1);
        gft_bl_true = w_hat_true_bl;
    elseif isfield(bl.p, 'w_hat_mse') && isfield(bl.p, 'w_hat_true')
        % Use pre-computed MSE if available
        gft_bl_mse = bl.p.w_hat_mse;
        gft_bl_true = bl.p.w_hat_true;
    else
        gft_bl_mse = [];
        gft_bl_true = [];
    end
    
    if isfield(qq.p, 'w_hat') && isfield(qq.p, 'w_hat_true')
        % Remove rows with non-zero imaginary parts in EITHER w_hat or w_hat_true
        w_hat_qq = qq.p.w_hat;
        w_hat_true_qq = qq.p.w_hat_true;
        
        % Filter out rows (countries) with non-zero imaginary parts in either variable
        tol = 1e-12;
        valid_qq = all(abs(imag(w_hat_qq)) <= tol, 2) & all(abs(imag(w_hat_true_qq)) <= tol, 2);
        w_hat_qq = real(w_hat_qq(valid_qq, :));
        w_hat_true_qq = real(w_hat_true_qq(valid_qq, :));
        
        % Check if we still have data after filtering
        if isempty(w_hat_qq) || sum(valid_qq) == 0
            warning('Spec %s (QQ): All countries have complex GFT values, skipping GFT plot', spec);
            gft_qq_mse = [];
            gft_qq_true = [];
        else
            % Compute MSE (mean over countries for each simulation)
            gft_qq_mse = mean((w_hat_qq - w_hat_true_qq).^2, 1);
            gft_qq_true = w_hat_true_qq;
        end
    elseif isfield(qq.p, 'w_hat_mse') && isfield(qq.p, 'w_hat_true')
        % Use pre-computed MSE if available
        gft_qq_mse = qq.p.w_hat_mse;
        gft_qq_true = qq.p.w_hat_true;
    else
        gft_qq_mse = [];
        gft_qq_true = [];
    end
    
    % Store data
    data = struct(...
        'lnn', lnn, ...
        'median_bl', median_bl, ...
        'median_qq', median_qq, ...
        'median_tr', median_tr, ...
        'bl_LB', bl_LB, ...
        'bl_UB', bl_UB, ...
        'qq_LB', qq_LB, ...
        'qq_UB', qq_UB, ...
        'gft_bl_mse', gft_bl_mse, ...
        'gft_bl_true', gft_bl_true, ...
        'gft_qq_mse', gft_qq_mse, ...
        'gft_qq_true', gft_qq_true);
end

function gmm_data = load_gmm_estimators(baseline_path, version, spec, estimators)
    % Load additional GMM estimators for comparison
    n_est = length(estimators);
    gmm_data = cell(n_est, 1);
    
    for i = 1:n_est
        est = estimators{i};
        basis = est.basis;
        knots = est.knots;
        
        % Construct filename
        gmm_file = fullfile(baseline_path, ...
            sprintf("Simulated_%s_Estimates_%s_%d_%s.mat", spec, basis, knots, version));
        
        if exist(gmm_file, 'file')
            gmm = load(gmm_file);
            
            % Compute support
            lnn = linspace(prctile(gmm.Rn_mat, 0.5, 'all'), ...
                           prctile(gmm.Rn_mat, 99.5, 'all'), 500)';
            
            % Extract medians and bands
            low = 2.5; high = 97.5;
            gmm_data{i} = struct(...
                'lnn', lnn, ...
                'median', median(gmm.p.theta_n, 2), ...
                'LB', prctile(gmm.p.theta_n, low, 2), ...
                'UB', prctile(gmm.p.theta_n, high, 2), ...
                'label', sprintf('%s (k=%d)', basis, knots), ...
                'basis', basis, ...
                'knots', knots);
        else
            warning('GMM file not found: %s', gmm_file);
            gmm_data{i} = [];
        end
    end
    
    % Remove empty cells
    gmm_data = gmm_data(~cellfun(@isempty, gmm_data));
end

function plot_elasticity_estimates(data, specs, spec_map, output_path, ...
                                   save_separate, save_combined, gmm_estimators, ylim_theta, ...
                                   file_suffix, legend_labels, legend_position, font_size, legend_font_size)
    n_specs = length(specs);
    
    % Set default text interpreter to LaTeX for this figure
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    
    % Align all specifications to common support
    data = align_to_common_support(data, specs);
    
    % Color palette
    pal.qq   = [0.60 0.00 0.00];
    pal.gmm1 = [0.00 0.20 0.55];      % Baseline GMM
    pal.gmm2 = [0.85 0.33 0.10];      % Additional GMM 1
    pal.gmm3 = [0.47 0.67 0.19];      % Additional GMM 2
    pal.true = [0.00 0.00 0.00];
    
    gmm_colors = {pal.gmm1, pal.gmm2, pal.gmm3};
    
    % Line styles
    sty.qq   = {'--', 'LineWidth', 3.5};
    sty.gmm1 = {':', 'LineWidth', 3.5};
    sty.gmm2 = {'-.', 'LineWidth', 3.0};
    sty.gmm3 = {'-', 'LineWidth', 3.5};  % Thicker for better visibility
    sty.true = {'-', 'LineWidth', 2.2};
    
    gmm_styles = {sty.gmm1, sty.gmm2, sty.gmm3};
    
    % X-axis settings
    xt = log(0.0001):log(10):0;
    xtlabs = arrayfun(@(x) sprintf('%.3g', exp(x)), xt, 'UniformOutput', false);
    
    % Layout parameters
    top_margin = 0.07; bottom = 0.12; left = 0.06; right = 0.04;
    gap = 0.05;
    width = (1 - left - right - (n_specs-1)*gap) / n_specs;
    height = 1 - bottom - top_margin;
    
    % Create figure
    f = figure('Color', 'w', 'Position', [100 100 300*n_specs+500 600]);
    
    axes_handles = gobjects(n_specs, 1);
    
    for i = 1:n_specs
        spec = specs{i};
        spec_data = data.(spec);
        
        % Create subplot
        ax = subplot(1, n_specs, i, 'Parent', f);
        axes_handles(i) = ax;
        hold(ax, 'on');
        
        % Set position
        pos_left = left + (i-1)*(width + gap);
        set(ax, 'PositionConstraint', 'outerposition', ...
                'OuterPosition', [pos_left, bottom, width, height]);
        
        % Plot baseline GMM confidence band
        shaded_ci(ax, spec_data.lnn_common, spec_data.bl_LB, spec_data.bl_UB, ...
                  gmm_colors{1}, 0.20);
        
        % Plot QQ confidence band
        shaded_ci(ax, spec_data.lnn_common, spec_data.qq_LB, spec_data.qq_UB, ...
                  pal.qq, 0.20);
        
        % Plot baseline GMM line
        h_gmm(1) = plot(spec_data.lnn_common, spec_data.median_bl, gmm_styles{1}{:}, ...
                    'Color', gmm_colors{1});
        
        % Determine labels - use custom if provided, otherwise use defaults
        if ~isempty(legend_labels)
            % User provided custom labels
            gmm_labels_to_use = legend_labels;
        else
            % Use default labels based on number of estimators
            if isfield(spec_data, 'gmm_extra') && ~isempty(spec_data.gmm_extra)
                % Multiple GMM estimators
                gmm_labels_to_use = {...
                    'Semiparametric GMM', ...
                    'Pareto GMM Estimator', ...
                    'Log-Normal GMM', ...
                    'Log-Normal QQ Estimator', ...
                    'True Function'};
            else
                % Single GMM estimator
                gmm_labels_to_use = {...
                    'Semiparametric GMM', ...
                    'Log-Normal QQ Estimator', ...
                    'True Function'};
            end
        end
        
        % Plot additional GMM estimators (if any)
        if isfield(spec_data, 'gmm_extra') && ~isempty(spec_data.gmm_extra)
            for j = 1:length(spec_data.gmm_extra)
                gmm_j = spec_data.gmm_extra{j};
                color_idx = min(j+1, length(gmm_colors));
                style_idx = min(j+1, length(gmm_styles));
                
                % Confidence band with lighter alpha
                shaded_ci(ax, gmm_j.lnn_common, gmm_j.LB, gmm_j.UB, ...
                          gmm_colors{color_idx}, 0.15);
                
                % Line
                h_gmm(j+1) = plot(gmm_j.lnn_common, gmm_j.median, ...
                                  gmm_styles{style_idx}{:}, 'Color', gmm_colors{color_idx});
            end
        end
        
        % Plot QQ line
        h_qq = plot(spec_data.lnn_common, spec_data.median_qq, sty.qq{:}, ...
                    'Color', pal.qq);
        
        % Plot true function
        h_tr = plot(spec_data.lnn_common, spec_data.median_tr, sty.true{:}, ...
                    'Color', pal.true);
        
        % Add legend only to first panel
        if i == 1
            legend_handles = [h_gmm, h_qq, h_tr];
            lgd = legend(legend_handles, gmm_labels_to_use, ...
                'Location', legend_position, 'Box', 'off','FontSize', legend_font_size, 'Interpreter', 'latex');
            set(lgd, 'FontName', 'Times');
        end
        
        % Formatting
        xticks(xt);
        xticklabels(xtlabs);
        xlabel('Log of Exporter Firm Share', 'FontWeight', 'normal', 'Interpreter', 'latex');
        ylabel('$\theta(\sigma-1)$', 'Interpreter', 'latex');
        ylim(ylim_theta);
        % title(spec_map.(spec).label, 'FontWeight', 'normal');

        ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];
        ax.XTickLabel = {'0.01\%' ,'0.1\%','1\%','10\%','100\%'};
        xlim([-9.5 -2]);
        % Set tick label interpreter to latex AFTER setting labels
        ax.TickLabelInterpreter = 'latex';

%         xlim([log(.0001) log(.1)])
        ax.FontSize = font_size;

        hold off;
    end
    
    % Apply formatting to all axes
    set(axes_handles, 'PositionConstraint', 'outerposition');
    set(axes_handles, 'LabelFontSizeMultiplier', 1.25, ...
                      'TitleFontSizeMultiplier', 1.25);
    tick_scale = 1.3;
    for ax = axes_handles'
        ax.FontSize = ax.FontSize * tick_scale;
    end
    
    % Overall title
    sgt = sgtitle(f, 'Trade Elasticity Function', 'FontWeight', 'bold', 'Interpreter', 'latex');
    sgt.FontSize = sgt.FontSize * 1.5;
    
    set(f, 'Renderer', 'painters');
    
    % Save files
    if save_combined
        if ~isempty(file_suffix)
            save_name = sprintf("Elasticity_%s_%s.pdf", strjoin(specs, '_'), file_suffix);
        else
            save_name = sprintf("Elasticity_%s.pdf", strjoin(specs, '_'));
        end
        print(f, fullfile(output_path, save_name), '-dpdf', '-painters');
    end
    
    if save_separate
        for i = 1:n_specs
            spec = specs{i};
            if ~isempty(file_suffix)
                pdf_out = fullfile(output_path, sprintf("Elasticity_%s_%s.pdf", spec, file_suffix));
            else
                pdf_out = fullfile(output_path, sprintf("Elasticity_%s.pdf", spec));
            end
            save_panel(axes_handles(i), pdf_out, i == 1, ylim_theta);
        end
    end
    
    close(f);
    
    % Reset default interpreters
    set(groot, 'defaultTextInterpreter', 'remove');
    set(groot, 'defaultLegendInterpreter', 'remove');
    set(groot, 'defaultAxesTickLabelInterpreter', 'remove');
end

function plot_gft_histograms(data, specs, spec_map, output_path, save_separate, xlim_gft, gft_tick_delta, gft_binwidth, font_size, legend_font_size)
    n_specs = length(specs);
    
    % Set default text interpreter to LaTeX for this figure
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    
    % Check if any spec has GFT data
    has_gft = false;
    for i = 1:n_specs
        spec = specs{i};
        if ~isempty(data.(spec).gft_bl_mse)
            has_gft = true;
            break;
        end
    end
    
    if ~has_gft
        warning('No GFT data found in loaded files. Skipping GFT plots.');
        return;
    end
    
    % Color palette
    pal.qq = [0.60 0.00 0.00];   % QQ estimator
    pal.base = [0.00 0.20 0.55]; % Baseline GMM
    
    % X-axis settings based on xlim_gft and gft_tick_delta
    % Calculate number of ticks from the range and tick spacing
    xt = xlim_gft(1):gft_tick_delta:xlim_gft(2);
    xtlabs = arrayfun(@(x) sprintf('%.3g', x), xt, 'UniformOutput', false);
    
    % Layout parameters
    top_margin = 0.07; bottom = 0.12; left = 0.06; right = 0.04;
    gap = 0.05;
    width = (1 - left - right - (n_specs-1)*gap) / n_specs;
    height = 1 - bottom - top_margin;
    
    % Create figure
    f = figure('Color', 'w', 'Position', [100 100 300*n_specs+500 600]);
    
    axes_handles = gobjects(n_specs, 1);
    
    for i = 1:n_specs
        spec = specs{i};
        spec_data = data.(spec);
        
        % Skip if no GFT data for this spec
        if isempty(spec_data.gft_bl_mse)
            continue;
        end
        
        % Create subplot
        ax = subplot(1, n_specs, i);
        axes_handles(i) = ax;
        hold on;
        
        % Set position
        pos_left = left + (i-1)*(width + gap);
        set(ax, 'PositionConstraint', 'outerposition', ...
                'OuterPosition', [pos_left, bottom, width, height]);
        
        % Compute RMSE / Average true GFT
        rmse_qq = sqrt(spec_data.gft_qq_mse) ./ mean(spec_data.gft_qq_true, 1);
        rmse_bl = sqrt(spec_data.gft_bl_mse) ./ mean(spec_data.gft_bl_true, 1);
        
        % % Debug output
        % fprintf('Spec %s: rmse_qq = %s, rmse_bl = %s\n', spec, mat2str(rmse_qq), mat2str(rmse_bl));
        % fprintf('  gft_qq_mse = %s, gft_bl_mse = %s\n', mat2str(spec_data.gft_qq_mse), mat2str(spec_data.gft_bl_mse));
        % fprintf('  mean(gft_qq_true) = %f, mean(gft_bl_true) = %f\n', ...
        %         mean(spec_data.gft_qq_true, 1), mean(spec_data.gft_bl_true, 1));
        
        fprintf('GFT Mean:\t%s: rmse_qq = %4.4f, rmse_bl = %4.4f\n', spec, mean(rmse_qq), mean(rmse_bl));
        fprintf('GFT Median:\t%s: rmse_qq = %4.4f, rmse_bl = %4.4f\n', spec, median(rmse_qq), median(rmse_bl));

        % Remove NaN or Inf values
        rmse_qq = rmse_qq(isfinite(rmse_qq));
        rmse_bl = rmse_bl(isfinite(rmse_bl));
        
        % Check if we have valid data
        has_qq = ~isempty(rmse_qq);
        has_bl = ~isempty(rmse_bl);
        
        if ~has_qq && ~has_bl
            continue;  % Skip this panel if no valid data
        end
        
        % Plot histograms
        if has_qq
            h_qq = histogram(rmse_qq, ...
                'FaceColor', pal.qq, ...
                'EdgeColor', 'none', ...
                'Normalization', 'probability', ...
                'FaceAlpha', 0.40, ...
                'BinWidth', gft_binwidth);
        end
        
        if has_bl
            h_bl = histogram(rmse_bl, ...
                'FaceColor', pal.base, ...
                'EdgeColor', 'none', ...
                'Normalization', 'probability', ...
                'FaceAlpha', 0.40, ...
                'BinWidth', gft_binwidth);
        end
        
        % Add legend only to first panel
        if i == 1
            if has_qq && has_bl
                lgd = legend([h_qq h_bl], ...
                    {'Log-normal QQ Estimator', 'Semiparametric GMM'}, ...
                    'Location', 'northeast', 'Box', 'off','FontSize', legend_font_size, 'Interpreter', 'latex');
                set(lgd, 'FontName', 'Times');
            elseif has_qq
                lgd = legend(h_qq, {'Log-normal QQ Estimator'}, ...
                    'Location', 'northeast', 'Box', 'off','FontSize', legend_font_size, 'Interpreter', 'latex');
                set(lgd, 'FontName', 'Times');
            elseif has_bl
                lgd = legend(h_bl, {'Semiparametric GMM'}, ...
                    'Location', 'northeast', 'Box', 'off','FontSize', legend_font_size, 'Interpreter', 'latex');
                set(lgd, 'FontName', 'Times');
            end
            if exist('lgd', 'var')
                % fontsize(lgd, 'scale', 1.5);
            end
        end
        
        % Formatting
        set(ax, 'XLim', xlim_gft, 'XTick', xt, 'XTickLabel', xtlabs);
        xlabel('MSE of Estimated GFT', 'FontWeight', 'normal', 'Interpreter', 'latex');
        ax.TickLabelInterpreter = 'latex';
        ax.FontSize = font_size;
        % title(spec_map.(spec).label, 'FontWeight', 'normal');
        
        hold off;
    end
    
    % Apply formatting
    set(axes_handles(isgraphics(axes_handles)), 'LabelFontSizeMultiplier', 1.25, ...
                      'TitleFontSizeMultiplier', 1.25);
    tick_scale = 1.3;
    for ax = axes_handles'
        if isgraphics(ax)
            ax.FontSize = ax.FontSize * tick_scale;
        end
    end
    
    % Overall title
    sgt = sgtitle(f, 'Gains from Trade: RMSE/Avg', 'FontWeight', 'bold', 'Interpreter', 'latex');
    sgt.FontSize = sgt.FontSize * 1.5;
    
    set(f, 'Renderer', 'painters');
    
    % Save files
    if save_separate
        for i = 1:n_specs
            if isgraphics(axes_handles(i))
                spec = specs{i};
                pdf_out = fullfile(output_path, ...
                    sprintf("GFT_%s.pdf", spec));
                save_panel(axes_handles(i), pdf_out, i == 1, []);
            end
        end
    end
    
    close(f);
    
    % Reset default interpreters
    set(groot, 'defaultTextInterpreter', 'remove');
    set(groot, 'defaultLegendInterpreter', 'remove');
    set(groot, 'defaultAxesTickLabelInterpreter', 'remove');
end

function data = align_to_common_support(data, specs)
    % Align all specifications to a common support using ismembertol
    % This ensures all panels are plotted on comparable domains
    % Also aligns GMM estimators within each spec to the common support
    
    if length(specs) == 1
        % Single specification - just rename lnn to lnn_common
        spec = specs{1};
        data.(spec).lnn_common = data.(spec).lnn;
        
        % Align extra GMM estimators to the main support
        if isfield(data.(spec), 'gmm_extra') && ~isempty(data.(spec).gmm_extra)
            for j = 1:length(data.(spec).gmm_extra)
                gmm_j = data.(spec).gmm_extra{j};
                [tf, idx_in_main] = ismembertol(gmm_j.lnn, data.(spec).lnn, 1e-2, 'DataScale', 1);
                idx_gmm_in_main = find(tf);
                idx_main_match = idx_in_main(tf);
                
                data.(spec).gmm_extra{j}.lnn_common = data.(spec).lnn(idx_main_match);
                data.(spec).gmm_extra{j}.median = gmm_j.median(idx_gmm_in_main);
                data.(spec).gmm_extra{j}.LB = gmm_j.LB(idx_gmm_in_main);
                data.(spec).gmm_extra{j}.UB = gmm_j.UB(idx_gmm_in_main);
            end
        end
        return;
    end
    
    tol = 1e-2;
    
    % Use first specification as reference
    ref_spec = specs{1};
    lnn_ref = data.(ref_spec).lnn;
    
    % Align each specification to the reference
    for i = 1:length(specs)
        spec = specs{i};
        lnn_curr = data.(spec).lnn;
        
        if i == 1
            % Reference specification - align to itself
            data.(spec).lnn_common = lnn_ref;
        else
            % Find common support between current and reference
            [tf, idx_in_ref] = ismembertol(lnn_curr, lnn_ref, tol, 'DataScale', 1);
            idx_curr_in_ref = find(tf);
            idx_ref_match = idx_in_ref(tf);
            
            % Store common support
            lnn_common = lnn_ref(idx_ref_match);
            
            % Filter all data to common support
            data.(spec).median_bl = data.(spec).median_bl(idx_curr_in_ref);
            data.(spec).median_qq = data.(spec).median_qq(idx_curr_in_ref);
            data.(spec).median_tr = data.(spec).median_tr(idx_curr_in_ref);
            data.(spec).bl_LB = data.(spec).bl_LB(idx_curr_in_ref);
            data.(spec).bl_UB = data.(spec).bl_UB(idx_curr_in_ref);
            data.(spec).qq_LB = data.(spec).qq_LB(idx_curr_in_ref);
            data.(spec).qq_UB = data.(spec).qq_UB(idx_curr_in_ref);
            data.(spec).lnn_common = lnn_common;
            
            % Also filter reference specification if first pass
            if i == 2
                data.(ref_spec).median_bl = data.(ref_spec).median_bl(idx_ref_match);
                data.(ref_spec).median_qq = data.(ref_spec).median_qq(idx_ref_match);
                data.(ref_spec).median_tr = data.(ref_spec).median_tr(idx_ref_match);
                data.(ref_spec).bl_LB = data.(ref_spec).bl_LB(idx_ref_match);
                data.(ref_spec).bl_UB = data.(ref_spec).bl_UB(idx_ref_match);
                data.(ref_spec).qq_LB = data.(ref_spec).qq_LB(idx_ref_match);
                data.(ref_spec).qq_UB = data.(ref_spec).qq_UB(idx_ref_match);
                data.(ref_spec).lnn_common = lnn_common;
            end
        end
        
        % Align extra GMM estimators within this spec to the common support
        if isfield(data.(spec), 'gmm_extra') && ~isempty(data.(spec).gmm_extra)
            for j = 1:length(data.(spec).gmm_extra)
                gmm_j = data.(spec).gmm_extra{j};
                [tf_gmm, idx_in_common] = ismembertol(gmm_j.lnn, data.(spec).lnn_common, tol, 'DataScale', 1);
                idx_gmm_in_common = find(tf_gmm);
                idx_common_match = idx_in_common(tf_gmm);
                
                data.(spec).gmm_extra{j}.lnn_common = data.(spec).lnn_common(idx_common_match);
                data.(spec).gmm_extra{j}.median = gmm_j.median(idx_gmm_in_common);
                data.(spec).gmm_extra{j}.LB = gmm_j.LB(idx_gmm_in_common);
                data.(spec).gmm_extra{j}.UB = gmm_j.UB(idx_gmm_in_common);
            end
        end
    end
end

function shaded_ci(ax, x, lb, ub, col, alpha)
    x = x(:); lb = lb(:); ub = ub(:);
    fill(ax, [x; flipud(x)], [lb; flipud(ub)], col, ...
         'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function save_panel(ax, pdf_out, keepLegend, ylim_theta)
    if nargin < 3, keepLegend = false; end
    if nargin < 4 || isempty(ylim_theta), ylim_theta = []; end
    
    ftmp = figure('Color', 'w', 'Position', [100 100 800 600]);
    newAx = copyobj(ax, ftmp);
    set(newAx, 'Units', 'normalized', 'Position', [0.13 0.11 0.775 0.815], ...
               'PositionConstraint', 'outerposition');
    
    % Set y-axis limits on the copied axis only if specified
    if ~isempty(ylim_theta)
        ylim(newAx, ylim_theta);
    end
    
    lgdsrc = legend(ax);
    if keepLegend && ~isempty(lgdsrc) && isvalid(lgdsrc)
        srcItems = lgdsrc.PlotChildren(:);
        strs = string(lgdsrc.String(:));
        tgtItems = gobjects(size(srcItems));
        for k = 1:numel(srcItems)
            s = srcItems(k);
            cands = findobj(newAx, 'Type', get(s, 'Type'));
            try
                props = {'LineStyle', 'LineWidth', 'Color', 'Marker', ...
                         'EdgeColor', 'FaceColor'};
                isMatch = false(size(cands));
                for j = 1:numel(cands)
                    same = true;
                    for p = 1:numel(props)
                        if isprop(s, props{p}) && isprop(cands(j), props{p})
                            vs = get(s, props{p});
                            vc = get(cands(j), props{p});
                            if ~isequal(vs, vc)
                                same = false;
                                break;
                            end
                        end
                    end
                    isMatch(j) = same;
                end
                idx = find(isMatch, 1, 'first');
                if ~isempty(idx), tgtItems(k) = cands(idx); end
            catch
                if ~isempty(cands), tgtItems(k) = cands(1); end
            end
        end
        good = isgraphics(tgtItems);
        if any(good)
            lgdnew = legend(newAx, tgtItems(good), cellstr(strs(good)), ...
                            'Location', lgdsrc.Location, 'Box', lgdsrc.Box);
            try, lgdnew.FontSize = lgdsrc.FontSize; end
        end
    end
    set(ftmp, 'Renderer', 'painters', 'Units', 'inches');
    figpos = get(ftmp, 'Position');
    marg = 0.02;
    set(ftmp, 'PaperUnits', 'inches', ...
        'PaperSize', [figpos(3)+2*marg, figpos(4)+2*marg], ...
        'PaperPosition', [marg, marg, figpos(3), figpos(4)], ...
        'InvertHardCopy', 'off');
    print(ftmp, pdf_out, '-dpdf', '-painters');
    close(ftmp);
end