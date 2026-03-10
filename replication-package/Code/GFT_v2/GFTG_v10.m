function [F,N_hat] = GFTG_v10(n_ii_hat,n_ij,x_ij,G_ij,denominator,GrhoA,GepsilonA,GratioA,Grho_barA,solveNi)


    x_ii = diag(x_ij);
    n_ii = diag(n_ij);
    G_ii = diag(G_ij);

    C = size(n_ii,1);

% eq:N_i_GT

if solveNi == 1
    gamma_part1 = arrayfun(@(n,G) GepsilonA(n,G)./n./GrhoA(n,G), n_ii.*n_ii_hat,G_ii);
    % bot_fint = @(lb,ub,G) integral(@(n) Grho_barA(n,G),lb,ub);
    % 
    % gamma_part1 = arrayfun(@(n,G) GepsilonA(n,G), n_ii.*n_ii_hat,G_ii) ...
    %     ./arrayfun(bot_fint,zeros(C,1),n_ii.*n_ii_hat,G_ii);


    top_fint = @(lb,ub,G) integral(@(n) GratioA(n,G),lb,ub); 
    N_hat = (1-gamma_part1.*arrayfun(top_fint,zeros(C,1),n_ii.*n_ii_hat,G_ii))./denominator;
else
    N_hat = ones(C,1);
end

F = x_ii.*n_ii_hat.*N_hat.*arrayfun(GrhoA,n_ii.*n_ii_hat,G_ii)./arrayfun(GrhoA,n_ii,G_ii) -arrayfun(GepsilonA,n_ii.*n_ii_hat,G_ii)./arrayfun(GepsilonA,n_ii,G_ii);

