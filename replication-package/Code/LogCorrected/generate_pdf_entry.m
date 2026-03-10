function [e_vals, pdf_vals] = generate_pdf_entry(e_bar, alpha_e, num_points)
    % Generates and plots the PDF g^e(e) for gamma^e = -1, 0, 1
    %
    % Inputs:
    %   e_bar: lower bound parameter (underline{e}), must be > 1
    %   alpha_e: alpha^e parameter, must be > 0
    %   num_points: number of points to evaluate (default: 1000)
    %
    % Outputs:
    %   e_vals: vector of e values
    %   pdf_vals: matrix of PDF values (rows = e values, columns = gamma values)
    %
    % Distribution: G^e(e) = 1 - (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e)
    % PDF: g^e(e) = (1/e) * (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e) 
    %              * [alpha_e + gamma_e/ln(e)]
    
    if nargin < 3
        num_points = 1000;
    end
    
    % Validate inputs
    if e_bar <= 1
        error('e_bar must be greater than 1');
    end
    if alpha_e <= 0
        error('alpha_e must be greater than 0');
    end
    
    % Generate e values over [e_bar, 10]
    e_vals = linspace(e_bar+.1, 5, num_points);
    
    % Gamma values to plot
    gamma_values = [-1, 0, 1];
    pdf_vals = zeros(num_points, length(gamma_values));
    
    % Create figure
    figure;
    hold on;
    
    % Calculate and plot PDF for each gamma value
    for i = 1:length(gamma_values)
        gamma_e = gamma_values(i);
        
        % Calculate PDF
        % g^e(e) = (1/e) * (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e) 
        %          * [alpha_e + gamma_e/ln(e)]
        
        term1 = 1 ./ e_vals;
        term2 = (e_vals / e_bar).^(-alpha_e);
        term3 = (log(e_vals) / log(e_bar)).^(-gamma_e);
        term4 = alpha_e + gamma_e ./ log(e_vals);
        
        pdf_vals(:, i) = term1 .* term2 .* term3 .* term4;
        
        % Plot
        plot(e_vals, pdf_vals(:, i), 'LineWidth', 2, 'DisplayName', ...
             sprintf('γ^e = %d', gamma_e));
        
        % Verify that PDF integrates to approximately 1
        integral_check = trapz(e_vals, pdf_vals(:, i));
        fprintf('Integral of PDF (γ^e = %d) over [%.2f, 10]: %.4f\n', ...
                gamma_e, e_bar, integral_check);
        if integral_check < 0.95
            warning('PDF integral for γ^e = %d is %.4f, which may indicate the range [e_bar, 10] does not capture most of the distribution', ...
                    gamma_e, integral_check);
        end
    end
    
    hold off;
    
%     % Add vertical line at e_bar
%     xline(e_bar, '--k', 'LineWidth', 1.5, 'Alpha', 0.5, ...
%           'DisplayName', sprintf('ē = %.2f', e_bar));
        ax = gca;
%     log_n_ticks = log([.0001 .001 0.01, 0.1, 1]);
%     ax.XTick = log_n_ticks;
%     ax.XTickLabel = {'.01%', '.1%','1%', '10%', '100%'};
    ax.FontSize = 16; 

    xlabel('e', 'FontSize', 20);
    ylabel('g^e(e)', 'FontSize', 20);
%     title(sprintf('PDF: ē = %.2f, α^e = %.2f', ...
%                   e_bar, alpha_e), 'FontSize', 18, 'Interpreter', 'tex');
    legend('Location', 'best', 'Interpreter', 'tex');
    grid on;
end