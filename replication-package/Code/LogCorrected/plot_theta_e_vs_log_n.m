function plot_thetae_vs_log_n(alpha_e, e_bar, n_min, n_max, num_points)
    % Plots theta^e(n) vs log(n) for gamma^e = -1, 0, 1
    %
    % Inputs:
    %   alpha_e: value of alpha^e parameter (must be > 0)
    %   e_bar: lower bound parameter (underline{e}), must be > 1
    %   n_min: minimum value of n (default: 0.01, must be in (0,1))
    %   n_max: maximum value of n (default: 0.99, must be in (0,1))
    %   num_points: number of points to evaluate (default: 100)
    %
    % Formula: theta^e(n) = alpha^e + gamma^e/ln(epsilon(n))
    % where epsilon(n) = (G^e)^{-1}(1-n)
    
    if nargin < 3
        n_min = 0.0001;
    end
    if nargin < 4
        n_max = 0.1;
    end
    if nargin < 5
        num_points = 1000;
    end
    
    % Validate inputs
    if e_bar <= 1
        error('e_bar must be greater than 1');
    end
    if alpha_e <= 0
        error('alpha_e must be greater than 0');
    end
    if n_min <= 0 || n_max >= 1
        error('n values must be in the interval (0, 1)');
    end
    
    % Generate n values
%     n = linspace(n_min, n_max, num_points);
        ln = linspace(log(n_min), log(n_max), num_points);
    n = exp(ln);

    % Calculate theta for different gamma values
    gamma_e_values = [-1, 0, 1];
    
    figure;
    hold on;
    
    for i = 1:length(gamma_e_values)
        gamma_e = gamma_e_values(i);
        theta_e = zeros(size(n));
        
        % For each n, compute epsilon(n) = (G^e)^{-1}(1-n)
        for j = 1:length(n)
            epsilon_n = inverse_cdf(1 - n(j), e_bar, alpha_e, gamma_e);
            theta_e(j) = alpha_e + gamma_e / log(epsilon_n);
        end
        
        plot(log(n), theta_e, 'LineWidth', 2, 'DisplayName', ...
             sprintf('γ^e = %d', gamma_e));
    end
    
    hold off;
    
    % Set custom x-axis labels with percentages
    ax = gca;
    log_n_ticks = log([.0001 .001 0.01, 0.1, 1]);
    ax.XTick = log_n_ticks;
    ax.XTickLabel = {'0.01%', '0.1%','1%', '10%', '100%'};
    ax.FontSize = 16; 
    xlabel('Log Exporter Firm Share ', 'FontSize', 20);
    ylabel('θ^e(n)', 'FontSize', 20, 'Interpreter', 'tex');
%     title(sprintf('θ^e(n) vs log(n) for α^e = %.2f, ē = %.2f', ...
%           alpha_e, e_bar), 'FontSize', 18, 'Interpreter', 'tex');
    legend('Location', 'northwest', 'Interpreter', 'tex', 'FontSize', 20);
    grid on;
    
%     % Add a horizontal line at alpha_e for reference
%     yline(alpha_e, '--k', 'LineWidth', 1, 'Alpha', 0.5, ...
%           'DisplayName', sprintf('α^e = %.2f', alpha_e));
end

function e = inverse_cdf(u, e_bar, alpha_e, gamma_e)
    % Solves G^e(e) = u for e
    % G^e(e) = 1 - (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e)
    
    % Define the equation to solve: G^e(e) - u = 0
    cdf_eq = @(e) 1 - (e/e_bar)^(-alpha_e) * (log(e)/log(e_bar))^(-gamma_e) - u;
    
    % Use fzero to find the root
    % Starting bounds: [e_bar, large value]
    try
        e = fzero(cdf_eq, [e_bar, 1e6]);
    catch
        % If fzero fails, try with a smaller upper bound
        e = fzero(cdf_eq, [e_bar, 1e3]);
    end
end