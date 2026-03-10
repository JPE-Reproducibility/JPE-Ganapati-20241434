function plot_theta_c_vs_log_n(alpha_c, n_min, n_max, num_points)
    % Plots theta^c(n) vs log(n) for gamma^c = -0.3, 0, 0.3
    %
    % Inputs:
    %   alpha_c: value of alpha^c parameter
    %   n_min: minimum value of n (default: 0.01, must be in (0,1))
    %   n_max: maximum value of n (default: 0.99, must be in (0,1))
    %   num_points: number of points to evaluate (default: 100)
    %
    % Formula: theta^c(n) = Er(n) / ((1/n) * integral_0^n Er(n') dn') - 1
    % where Er(n) = n^(-alpha^c) * (1 - ln(n))^(-gamma^c)
    
    if nargin < 2
        n_min = 0.0001;
    end
    if nargin < 3
        n_max = 0.1;
    end
    if nargin < 4
        num_points = 100;
    end
    
    % Validate inputs
    if n_min <= 0 || n_max >= 1
        error('n values must be in the interval (0, 1)');
    end
    
    % Generate n values
    ln = linspace(log(n_min), log(n_max), num_points);
    n = exp(ln);
%     n = linspace(n_min, n_max, num_points);
    
    % Gamma values to plot
    gamma_values = [-0.3, 0, 0.3];
    
    figure;
    hold on;
    
    for i = 1:length(gamma_values)
        gamma_c = gamma_values(i);
        theta_c = zeros(size(n));
        
        % For each n, compute theta^c(n)
        for j = 1:length(n)
            n_val = n(j);
            
            % Compute Er(n)
            Er_n = n_val^(-alpha_c) * (1 - log(n_val))^(-gamma_c);
            
            % Compute average: (1/n) * integral_0^n Er(n') dn'
            % Use numerical integration from small epsilon to n
            Er_func = @(x) x.^(-alpha_c) .* (1 - log(x)).^(-gamma_c);
            integral_val = integral(Er_func, 0, n_val, 'RelTol', 1e-6);
            avg_Er = integral_val / n_val;
            
            % Compute theta^c(n)
            theta_c(j) = Er_n / avg_Er - 1;
        end
        
        plot(log(n), theta_c, 'LineWidth', 2, 'DisplayName', ...
             sprintf('γ^c = %.1f', gamma_c));
    end
    
    hold off;
    
    % Set custom x-axis labels with percentages
    ax = gca;
    log_n_ticks = log([.0001 .001 0.01, 0.1, 1]);
    ax.XTick = log_n_ticks;
    ax.XTickLabel = {'0.01%', '0.1%','1%', '10%', '100%'};
    ax.FontSize = 16; 
    xlabel('Log Exporter Firm Share ', 'FontSize', 20);
    ylabel('θ^c(n)', 'FontSize', 20, 'Interpreter', 'tex');
%     title(sprintf('θ^c(n) vs log(n) for α^c = %.2f', alpha_c), ...
%           'FontSize', 18, 'Interpreter', 'tex');
    legend('Location', 'northwest', 'Interpreter', 'tex', 'FontSize', 20);
    grid on;
%     
%     % Add a horizontal line at 0 for reference
%     yline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5);
end