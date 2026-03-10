function plot_theta_combined(alpha_c, alpha_e, e_bar, n_min, n_max, num_points)
    % Plots theta(n) = 1 + (theta^c(n) + 1) * theta^e(n) vs log(n)
    % for gamma^c = -0.3, 0, 0.3
    %
    % Inputs:
    %   alpha_c: value of alpha^c parameter
    %   alpha_e: value of alpha^e parameter (must be > 0)
    %   e_bar: lower bound parameter (underline{e}), must be > 1
    %   n_min: minimum value of n (default: 0.01, must be in (0,1))
    %   n_max: maximum value of n (default: 0.99, must be in (0,1))
    %   num_points: number of points to evaluate (default: 100)
    %
    % Formula: theta(n) = 1 + (theta^c(n) + 1) * theta^e(n)
    % where:
    %   theta^c(n) = Er(n) / ((1/n) * integral_0^n Er(n') dn') - 1
    %   Er(n) = n^(-alpha^c) * (1 - ln(n))^(-gamma^c)
    %   theta^e(n) = alpha^e + gamma^e/ln(epsilon(n))
    %   epsilon(n) = (G^e)^{-1}(1-n)
    
    if nargin < 4
        n_min = 0.0001;
    end
    if nargin < 5
        n_max = 0.1;
    end
    if nargin < 6
        num_points = 100;
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
    ln = linspace(log(n_min), log(n_max), num_points);
    n = exp(ln);
%     n = linspace(n_min, n_max, num_points);
    
    % Gamma^c values to plot
    gamma_c_values = [-0.3, 0, 0.3];
    gamma_e_values = [-1, 0, 1];
    
    figure;
    hold on;
    
    for i = 1:length(gamma_c_values)
        gamma_c = gamma_c_values(i);
        gamma_e = gamma_e_values(i);

        theta_combined = zeros(size(n));
        
        % For each n, compute combined theta(n)
        for j = 1:length(n)
            n_val = n(j);
            
            % Compute theta^c(n)
            Er_n = n_val^(-alpha_c) * (1 - log(n_val))^(-gamma_c);
            Er_func = @(x) x.^(-alpha_c) .* (1 - log(x)).^(-gamma_c);
            integral_val = integral(Er_func, 0, n_val, 'RelTol', 1e-6);
            avg_Er = integral_val / n_val;
            theta_c = Er_n / avg_Er - 1;
            
            % Compute theta^e(n)
            epsilon_n = inverse_cdf(1 - n_val, e_bar, alpha_e, gamma_e);
            theta_e = alpha_e + gamma_e / log(epsilon_n);
            
            % Compute combined theta(n)
            theta_combined(j) = 1 + (theta_c + 1) * theta_e;
        end
        
        plot(log(n), theta_combined, 'LineWidth', 2, 'DisplayName', ...
             sprintf('γ^c = %.1f, γ^e = %.1f', gamma_c, gamma_e));
    end
    
    hold off;
    
    % Set custom x-axis labels with percentages
    ax = gca;
    log_n_ticks = log([.0001 .001 0.01, 0.1, 1]);
    ax.XTick = log_n_ticks;
    ax.XTickLabel = {'0.01%', '0.1%','1%', '10%', '100%'};
    ax.FontSize = 16; 
    xlabel('Log Exporter Firm Share ', 'FontSize', 20);
    ylabel('θ(n)', 'FontSize', 20, 'Interpreter', 'tex');
%     title(sprintf('θ(n) vs log(n): α^c = %.2f, α^e = %.2f, ē = %.2f', ...
%           alpha_c, alpha_e,  e_bar), 'FontSize', 18, 'Interpreter', 'tex');
    legend('Location', 'northwest', 'Interpreter', 'tex', 'FontSize', 20);
    grid on;
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