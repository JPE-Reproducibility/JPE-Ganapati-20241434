function [r2_vec,Ests_Baseline,B_vec,SE_vec] = QuantilesColombia(pctiles_raw,cty,plot)

N = sum(pctiles_raw(:,3));
n_ij_data = pctiles_raw(:,3)/N;
pctiles_data = pctiles_raw(cty,5:end);
n_ij        = n_ij_data(cty);

% Get this from the data for an ij pair (stanislav's code should have these
% values
delta_ij    = - 7           ;
nu_ij       = 1             ;% Variation in Sales

%% Baseline
% Optimize over delta and nu to fit the quantiles
p_grid      = .01:.01:.99;
F = @(x) QuantilesBaselineR2(x(1),x(2),n_ij,pctiles_data);
[Ests_Baseline]  = fminunc(@(x) abs(1-F(x))   ,[delta_ij,nu_ij],optimset('TolFun',1e-10,'TolX',1e-8,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
Quantiles_opt = QuantilesBaseline(Ests_Baseline(1),Ests_Baseline(2),n_ij);
Quantiles= Quantiles_opt;

% Get R2 for baseline
[B_baseline,SE_baseline,~,~,STATS] = regress(Quantiles,[ones(size(Quantiles)) pctiles_data']);
r2_Baseline = STATS(1);
Quantiles_pred_base = B_baseline(1)+pctiles_data*B_baseline(2);


if plot == 1
    figure;
    
    subplot(2,2,1); 
    scatter(Quantiles,p_grid)
    title('Baseline Cdf')
    
    subplot(2,2,2); 
    scatter(pctiles_data,p_grid)
    title('Empirical CDF')
    
    subplot(2,1,2); 
    scatter(pctiles_data,(Quantiles))
    title('Baseline vs Data')
    hold on
    line(pctiles_data,Quantiles_pred_base, 'Color', 'red')
    hold off
    saveas(gcf,'../../Output/QQ_Colombia_Baseline','epsc')

    sgtitle({'Baseline Economy fitted to Colombia-USA Trade',' n_{ij} = exporters to USA/total Exporters'})
    saveas(gcf,'../../Output/QQ_Colombia_Baseline','pdf')
end



%% Log-Normal
% Do quantiles under n = 1
m_ij_LN    = @(n) norminv(1 - n*((1 - p_grid))); %testing
n_ij_LN = 1;
m_ij_LN(n_ij_LN);

[B_LN,SE_LN,~,~,STATS] = regress(m_ij_LN(n_ij_LN)',[ones(size(Quantiles)) pctiles_data']);
r2_LN = STATS(1);
Quantiles_pred_LN = B_LN(1)+pctiles_data*B_LN(2);

if plot == 1
    
    figure;
    sgtitle({'Log Normal Economy Fitted to Colombia-US Trade',' n_{ij} = 1'})
    
    subplot(2,2,1); 
    scatter(m_ij_LN(n_ij_LN)',p_grid)
    title('lognormal m_{ij}')
    
    subplot(2,2,2); 
    scatter(pctiles_data,p_grid)
    title('Empirical CDF')
    
    subplot(2,1,2); 
    scatter(pctiles_data,m_ij_LN(n_ij_LN)')
    title('LogNormal vs Data')
    hold on
    line(pctiles_data,Quantiles_pred_LN, 'Color', 'red')
    hold off
    saveas(gcf,'../../Output/QQ_Colombia_LogNormal','pdf')

end

%% Log-Normal Under Selection
% Do quantiles under n = n_ij 
m_ij_LNS    = @(n) norminv(1 - n*((1 - p_grid))); %testing
n_ij_LNS = n_ij;
m_ij_LNS(n_ij_LNS);

[B_LNS,SE_LNS,~,~,STATS] = regress(m_ij_LNS(n_ij_LNS)',[ones(size(Quantiles)) pctiles_data']);
r2_LNS = STATS(1);
Quantiles_pred_LNS = B_LNS(1)+pctiles_data*B_LNS(2);

if plot == 1
    
    figure;
    sgtitle({'Log Normal Economy Fitted to Colombia-US Trade',' n_{ij} = data'})
    
    subplot(2,2,1); 
    scatter(m_ij_LNS(n_ij_LNS)',p_grid)
    title('lognormal m_{ij}')
    
    subplot(2,2,2); 
    scatter(pctiles_data,p_grid)
    title('Empirical CDF')
 
    subplot(2,1,2); 
    scatter(pctiles_data,m_ij_LNS(n_ij_LNS)')
    title('LogNormal vs Data')
    hold on
    line(pctiles_data,Quantiles_pred_LNS, 'Color', 'red')
    hold off
    saveas(gcf,'../../Output/QQ_Colombia_LogNormal_Selection','pdf')

end

%% Pareto
m_ij_Pareto = @(n) -log(1 - p_grid) + 0*n; 
n_ij_Pareto = 1;
m_ij_Pareto(n_ij_Pareto)

[B_pareto,SE_pareto,~,~,STATS] = regress(m_ij_Pareto(n_ij_Pareto)',[ones(size(Quantiles)) pctiles_data']);
r2_Pareto = STATS(1);
Quantiles_pred_pareto = B_pareto(1)+pctiles_data*B_pareto(2);

if plot == 1
    figure;
    sgtitle({'Pareto Economy Fitted to Colombia-US Trade',' n_{ij} = 1'})
    
    subplot(2,2,1); 
    scatter(m_ij_Pareto(n_ij_Pareto)',p_grid)
    title('Pareto m_{ij}')
    
    subplot(2,2,2); 
    scatter(pctiles_data,p_grid)
    title('Empirical CDF')
    
    subplot(2,1,2); 
    scatter(pctiles_data,m_ij_Pareto(n_ij_Pareto)')
    title('Pareto vs Data')
    hold on
    line(pctiles_data,Quantiles_pred_pareto, 'Color', 'red')
    hold off
    saveas(gcf,'../../Output/QQ_Colombia_Pareto','pdf')

end

if plot == 1
    figure;
    subplot(2,2,1); 
        scatter(pctiles_data,(Quantiles))
        title('Baseline vs Data')
        hold on
        line(pctiles_data,Quantiles_pred_base, 'Color', 'red')
        hold off
        set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        xlabel('Data Quantile');ylabel('Fitted Quantile');
    subplot(2,2,2); 
        scatter(pctiles_data,m_ij_LN(n_ij_LN)')
        title('LogNormal vs Data')
        hold on
        line(pctiles_data,Quantiles_pred_LN, 'Color', 'red')
        hold off
        set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        xlabel('Data Quantile');ylabel('Fitted Quantile');
    subplot(2,2,3); 
        scatter(pctiles_data,m_ij_LNS(n_ij_LNS)')
        title('LogNormal (Selection) vs Data')
        hold on
        line(pctiles_data,Quantiles_pred_LNS, 'Color', 'red')
        hold off
        set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        xlabel('Data Quantile');ylabel('Fitted Quantile');
    subplot(2,2,4); 
        scatter(pctiles_data,m_ij_Pareto(n_ij_Pareto)')
        title('Pareto vs Data')
        hold on
        line(pctiles_data,Quantiles_pred_pareto, 'Color', 'red')
        hold off
        set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        xlabel('Data Quantile');ylabel('Fitted Quantile');
            saveas(gcf,'../../Output/QQ_Colombia_Quantiles','pdf')

end

    figure;
        scatter(pctiles_data,(Quantiles))
%         title('Baseline vs Data')
        hold on
        line(pctiles_data,Quantiles_pred_base, 'Color', 'red')
        hold off
        set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        xlabel('Data Quantile','FontSize',18);ylabel('Fitted Quantile','FontSize',18);
        saveas(gcf,'../../Output/QQ_Colombia_Quantiles_BaselineP','epsc')

    figure;
        scatter(pctiles_data,m_ij_LN(n_ij_LN)')
%         title('LogNormal vs Data')
        hold on
        line(pctiles_data,Quantiles_pred_LN, 'Color', 'red')
        hold off
        set(gca,'YTickLabel',[]);set(gca,'XTickLabel',[]);
        xlabel('Data Quantile','FontSize',18);ylabel('Fitted Quantile','FontSize',18);
        saveas(gcf,'../../Output/QQ_Colombia_Quantiles_LNP','epsc')



    r2_vec = [r2_Baseline r2_LN r2_Pareto r2_LNS];
    B_vec  = [B_baseline(2) B_LN(2) B_pareto(2) B_LNS(2)];
    SE_vec = [abs(SE_baseline(2,1)-B_baseline(2))/1.96 abs(SE_LN(2,1)-B_LN(2))/1.96 abs(SE_pareto(2,1)-B_pareto(2))/1.96 abs(SE_LNS(2,1)-B_LNS(2))/1.96];
end