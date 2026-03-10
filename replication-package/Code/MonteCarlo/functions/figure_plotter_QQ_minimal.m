function figure_plotter_QQ_minimal(output_path, l, ii, jj, sim_id, fct_form, basis, knots)
    % plots predicted quantiles against the true ones from QQ regression
    % for a specific i-j pair
    
    % Extract data from l struct
    if ~isfield(l, 'y_within') || ~isfield(l, 'X1') || ~isfield(l, 'c') || ~isfield(l, 'pp')
        error('l struct must contain y_within, X1, c, and pp fields');
    end
    
    c = l.c;
    pp = l.pp;
    
    if numel(l.y_within) ~= c^2 * pp
        error('Q_lnx must be [pp*c^2 x 1].');
    end
    
    pair_idx = sub2ind([c c], ii, jj);
    Qmat = reshape(l.y_within, [c^2, pp]);   % rows: OD, cols: quantiles
    
    % extract data for this specific pair
    Q_true = Qmat(pair_idx, :).';               % [pp x 1]
    fitted = l.beta0 + l.beta * l.X1(:, 2);
    Qhatmat = reshape(fitted, [c^2, pp]);
    Q_pred = Qhatmat(pair_idx, :).';
    
    % set plots
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    
    % figure - match figures_pres_alt aspect ratio (800x600 for single panel)
    f = figure('Color', 'w', 'Position', [100 100 800 600]);
    ax = axes('Parent', f);
    
    % Set position to match figures_pres_alt single panel layout
    set(ax, 'Units', 'normalized', 'Position', [0.13 0.11 0.775 0.815], ...
            'PositionConstraint', 'outerposition');
    
    hold(ax, 'on');
    
    % color palette (matching figures_pres_alt)
    scatter_color = [0.00 0.20 0.55];  % Blue
    line_color = [0.60 0.00 0.00];     % Red
    
    % plot scatter
    scatter(Q_pred, Q_true, 36, scatter_color, 'filled', ...
            'MarkerFaceAlpha', 0.6);
    
    % add 45-degree reference line
    xlims = xlim(ax);
    ylims = ylim(ax);
    plot_lims = [min(xlims(1), ylims(1)), max(xlims(2), ylims(2))];
    plot(plot_lims, plot_lims, '-', 'Color', line_color, 'LineWidth', 2.2);
    
    % axes limits
    xlim(plot_lims);
    ylim(plot_lims);
    
    % formatting
    xlabel('Fitted Firm Sales', 'FontWeight', 'normal', 'Interpreter', 'latex');
    ylabel('Simulated Firm Sales', 'FontWeight', 'normal', 'Interpreter', 'latex');
    
    % font sizes matching figures_pres_alt style
    ax.FontSize = 12;
    ax.FontSize = ax.FontSize * 1.3;  % Apply tick scale
    
    % grid and box
    % grid off;
    % box on;
    
    % Title with R2 (optional)
    if isfield(l, 'R2_within')
        text(ax, 0.9, 0.1, sprintf('$R^2 = %.3f$', l.R2_within), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', ...
            'Interpreter', 'latex', ...
            'FontSize', ax.FontSize);
    end
    
    set(f, 'Renderer', 'painters');
    hold off;
    
    % save with proper paper size to match figure aspect ratio
    if nargin < 7 || isempty(knots)
        knots = 1;
    end
    
    % Set paper size to match figure size for PDF output
    set(f, 'Units', 'inches');
    figpos = get(f, 'Position');
    marg = 0.02;
    set(f, 'PaperUnits', 'inches', ...
        'PaperSize', [figpos(3)+2*marg, figpos(4)+2*marg], ...
        'PaperPosition', [marg, marg, figpos(3), figpos(4)], ...
        'InvertHardCopy', 'off');
    
    save_name = sprintf('QQ_scatter_%s_%s_k%d_sim%d_i%d_j%d.pdf', ...
                        fct_form, basis, knots, sim_id, ii, jj);
    pdf_out = fullfile(output_path, save_name);
    print(f, pdf_out, '-dpdf', '-painters');
    close(f);
    
    % reset default interpreters
    set(groot, 'defaultTextInterpreter', 'remove');
    set(groot, 'defaultLegendInterpreter', 'remove');
    set(groot, 'defaultAxesTickLabelInterpreter', 'remove');
end