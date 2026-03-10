function [Quantiles] = QuantilesBaseline(delta_ij,nu_ij,n_ij)
% Compute quantiles using baseline model

load baseline_estimate.mat
lnn     = @(n) evknots_make(o.k,log(n)) ;
lnd     = @(n) evknots_make_deriv(o.k,log(n)) ;
rho_bar = @(n) exp((lnn(n)*o.est_rho')') .* ((lnd(n)*o.est_rho')'+1);

% Get this from the data for an ij pair (stanislav's code should have these
% values
% n_ij        = n_ij_data(cty);
% delta_ij    = - 7;
% nu_ij       = 1   ;% Variation in Sales

% W need a distribution of sales that captures the support
lnx         = -6:.01:4.0;
p_grid      = .01:.01:.99;

% Define functions for CDF calculation
innerfct = @(x,n,n_ij,delta_ij,nu_ij) (x - log(rho_bar(n))')/nu_ij + delta_ij;
innerfct1 = @(x,n) normcdf(innerfct(x,n,n_ij,delta_ij,nu_ij))';
lnH = @(x,n_ij)  arrayfun(@(x,l,u) 1./u .* integral( @(n) innerfct1(x,n) ,l,u),x,zeros(size(n_ij)),n_ij);
CDF = lnH(exp(lnx),n_ij*ones(size(lnx)));

% Find unique values in CDF and their corresponding indices
[unique_cdf, ia, ~] = unique(CDF);

% Use these indices to get the corresponding unique X values
unique_lnx = lnx(ia);

% Get the quantiles
Quantiles = interp1(unique_cdf, unique_lnx, p_grid, 'spline', 'extrap')';

end