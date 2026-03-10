function plot_Er_vs_log_n(alpha_c, n_min, n_max, num_points)
    % Plots Er(n) vs log(n) for gamma^c = -0.3, 0, 0.3
    %
    % Inputs:
    %   alpha_c: value of alpha^c parameter
    %   n_min: minimum value of n (default: 0.01, must be in (0,1))
    %   n_max: maximum value of n (default: 0.99, must be in (0,1))
    %   num_points: number of points to evaluate (default: 500)
    %
    % Formula: Er(n) = n^(-alpha^c) * (1 - ln(n))^(-gamma^c)
    
    if nargin < 2
        n_min = 0.0001;
    end
    if nargin < 3
        n_max = 0.99;
    end
    if nargin < 4
        num_points = 500;
    end
    
    % Validate inputs
    if n_min <= 0 || n_max >= 1
        error('n values must be in the interval (0, 1)');
    end
    
    % Generate n values
    n = exp(linspace(log(n_min), log(n_max), num_points));
    
    % Gamma values to plot
    gamma_values = [-0.3, 0, 0.3];
    
    figure;
    hold on;
    
    for i = 1:length(gamma_values)
        gamma_c = gamma_values(i);
        
        % Compute Er(n) = n^(-alpha_c) * (1 - ln(n))^(-gamma_c)
        Er_n = n.^(-alpha_c) .* (1 - log(n)).^(-gamma_c);
        
        plot(log(n), Er_n, 'LineWidth', 2, 'DisplayName', ...
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
    ylabel('E[r|\epsilon(n)]', 'FontSize', 20, 'Interpreter', 'tex');
%     title(sprintf('E[r(n|e)] vs log(n) for α^c = %.2f', alpha_c), ...
%           'FontSize', 18, 'Interpreter', 'tex');
    legend('Location', 'best', 'Interpreter', 'tex');
    grid on;
end