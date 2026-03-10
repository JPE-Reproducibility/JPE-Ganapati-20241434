function [F,N_hat] = GFT(n_ii_hat,X_ij,epsilon,rho,rho_bar,ratio,n_ij,x_ij,kappa_i,precalculated_integrals,array,N_i_A_hat_flag,buffer)
%GFT Computes the Gains from Trade
%   There are two options
%    array -> if equals 'array' compute results simultaneiously for an full
%             list of countries
%             "The option ?array? computes GFT for an array of countries."
%    N_i_A_hat_flag -> if equals one compute results allowing N_i to vary,
%                      if not equal to 1, then fixes N_i
%                      "If the last entry N_i_A_hat_flag is not 1, then we 
%                      use the Pareto prediction that hatN = 1, 
%                      and modify the objective function."

    x_ii = diag(x_ij);
    n_ii = diag(n_ij);
    
    % this is the slowest step - approximations can help!
    % closed forms for pareto
    %,'RelTol',0,'AbsTol',1e-12
    top_fint = @(lb,ub) integral(ratio,lb,ub); 
    bot_fint = @(lb,ub) integral(rho,lb,ub);   

    switch array
        case 'array'
            top2 = arrayfun(top_fint,zeros(size(n_ii))+buffer,n_ii.*n_ii_hat);
            bot2 = arrayfun(bot_fint,zeros(size(n_ii))+buffer,n_ii.*n_ii_hat)./arrayfun(epsilon,n_ii.*n_ii_hat);
            N_i_A_hat = @(n_ii_hat) arrayfun(rho_bar,n_ii)./arrayfun(epsilon,n_ii).*arrayfun(epsilon,n_ii.*n_ii_hat)./arrayfun(rho_bar,n_ii.*n_ii_hat)./x_ii./n_ii_hat.*kappa_i;

            if N_i_A_hat_flag == 1
                F = sum(X_ij.*(1-precalculated_integrals),2)./sum(X_ij,2)-1./N_i_A_hat(n_ii_hat).*(1-top2./bot2);
            else
                F = 1 - N_i_A_hat(n_ii_hat);
            end
        otherwise
            % pre-compute these integrals
            top2 = arrayfun(top_fint,zeros(size(n_ii))+buffer,n_ii.*n_ii_hat);
            bot2 = arrayfun(bot_fint,zeros(size(n_ii))+buffer,n_ii.*n_ii_hat)./arrayfun(epsilon,n_ii.*n_ii_hat);
            N_i_A_hat = @(n_ii_hat) rho_bar(n_ii)./epsilon(n_ii).*epsilon(n_ii.*n_ii_hat)./rho_bar(n_ii.*n_ii_hat)./x_ii./n_ii_hat.*kappa_i;

            if N_i_A_hat_flag == 1
                F = sum(X_ij.*(1-precalculated_integrals),2)./sum(X_ij,2)-1./N_i_A_hat(n_ii_hat).*(1-top2./bot2);
            else
                F = 1- N_i_A_hat(n_ii_hat);
            end
    end
    
    N_hat = N_i_A_hat(n_ii_hat);

end

