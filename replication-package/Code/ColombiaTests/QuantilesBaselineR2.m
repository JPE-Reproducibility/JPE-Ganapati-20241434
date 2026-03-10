function [r2_Baseline] = QuantilesBaselineR2(delta_ij,nu_ij,n_ij,pctiles_data)
% Compute quantiles using baseline model and calculate R2

    try
        Quantiles = QuantilesBaseline(delta_ij,nu_ij,n_ij);
        [~,~,~,~,STATS] = regress(Quantiles,[ones(size(Quantiles)) pctiles_data']);
        r2_Baseline = STATS(1);
    catch
        r2_Baseline = 0;
    end

end