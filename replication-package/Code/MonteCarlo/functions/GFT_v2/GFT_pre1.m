function [F,N_hat] = GFT_pre1(index,n_ii_hat,X_ij,epsilon,rho,rho_bar,ratio,n_ij,x_ij,precalculated_integrals,top_fint,bot_fint)
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
    

    top2 = top_fint(0,n_ii(index).*n_ii_hat);
    bot2 = bot_fint(0,n_ii(index).*n_ii_hat)./epsilon(n_ii(index).*n_ii_hat);

    N_i_A_hat = @(n_ii_hat) rho_bar(n_ii(index)) ...
        ./epsilon(n_ii(index)).*epsilon(n_ii(index).*n_ii_hat) ...
        ./rho_bar(n_ii(index).*n_ii_hat)./x_ii(index)./n_ii_hat;
    
    F = sum(X_ij.*(1-precalculated_integrals),2)./sum(X_ij,2)-1./N_i_A_hat(n_ii_hat).*(1-top2./bot2);
    F = F(index);

    N_hat = N_i_A_hat(n_ii_hat);

end

