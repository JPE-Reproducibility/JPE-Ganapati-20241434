clear; close all; clc;
addpath('../GMM_estimation_v3')
load baseline_estimate.mat

lnn         = @(n) evknots_make(o.k,log(n)) ;
lnd         = @(n) evknots_make_deriv(o.k,log(n)) ;
epsilon     = @(n) exp((lnn(n)*o.est_epsilon')');
sigma       = 3.2;

ln_range    = -12:.01:-.01;
n_ij        = exp(ln_range);
elast_epsd  = diff(log(epsilon(n_ij)))./diff(ln_range);
N           = size(elast_epsd,2);

alpha       = .003;
beta        = 1.8;
gamma       = -1.0;
e_low       = -10;
e_high      = -4;

% F = @(x) match_LogPareto(x(1),x(2),x(3),x(4),x(5),elast_epsd,ln_range);
start = [e_high,e_low,alpha,beta,gamma];


% F2 = @(x) match_LogPareto(-4,-8.25,x(3),x(4),x(5),elast_epsd,ln_range);
F2 = @(x) match_LogPareto(-4,-8.25,x(3),x(4),x(5),elast_epsd,ln_range);
[sol,val] = fminunc(F2,start)
[FF2,elast_sol] = F2(sol);


plot(log((n_ij(1:N))),(1-sigma)*1./elast_sol,'-','LineWidth',3)
hold on
plot(log((n_ij(1:N))),(1-sigma)*1./elast_epsd,':','LineWidth',3)
legend('Baseline Estimate', 'Piecewise Log-Pareto','FontSize',18,'Location', 'northwest')
hold off
% ylabel('\partial ln \epsilon (n ) / \partial ln (n)','FontSize',18)
ylabel('$\theta^{e}(\sigma-1)$','FontSize',24,'Interpreter','Latex')

xlabel('Log Exporter Firm Share','FontSize',18)
ax = gca;
ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];
ax.XTickLabel = {'0.01%' ,'0.1%','1%','10%','100%'};
ax.FontSize = 18; 
box on
legend('boxoff')
saveas(gca,'../../Output/FA_logPareto_Sim','epsc')
