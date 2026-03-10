
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Pareto with log-normal entry (EKK-lite)   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % this only works with a uniform distribution (symmetric)
    xx = -9.3:.01:-.1;
    unif = exp(xx);
    unif_EKK = unif;
    sigma_lit       = 4;    %Bernard  et  al.  (2003).
    rng(1)
    r_vec           = zeros(size(unif_EKK));
    e_vec           = zeros(size(unif_EKK));
    iter            = 10000;
    theta_EKK       = 4;
    for i=1:1:iter
        EKK_mu      = [0 0];
        EKK_sigma   = [1.69 -.65; -.65 .34]; % Eaton, J., Kortum, S., and Kramarz, F. (2011). An anatomy of international trade: Evidence from french firms. Econometrica, 79(5):1453--1498
        R           = mvnrnd(EKK_mu,EKK_sigma,size(unif_EKK,2));
        z           = unif_EKK.^(-1/theta_EKK);
        EKK_B       = exp(R(:,1))';
        EKK_F       = exp(R(:,2))';
        r           = EKK_B.*(tau./z).^(1-sigma_lit);
        e           = r./EKK_F;
        [~,B]       = sort(e,'desc');
        e_vec       = e(B) + e_vec;
        r_vec       = r(B) + r_vec;
    end

    % Need to intialize so not infinity
    r_vec(1) = r_vec(2)*100000000;
    elast_r_EKK =     (r_vec/iter)./cumtrapz(unif_EKK,(r_vec/iter))./unif_EKK-1;
    e               = e_vec/iter;
    elast_e_EKK     = (log(e(1:(size(unif_EKK,2)-1)))-log(e(2:size(unif_EKK,2))))./(log(unif_EKK(1:(size(unif_EKK,2)-1)))-log(unif_EKK(2:size(unif_EKK,2))));
    r               = r_vec/iter;

    range = 1:1:900;
    range2 = 200:min(1400,size(unif_EKK,2)-1);
    ln_ij = log(unif(range));
    p=plot(log(unif_EKK(range2)),smoothdata(elast_e_EKK(range2),'gaussian',10),'-.'  )  ;

    r_smooth = smoothdata(r(range2)/iter,'gaussian',10)
    e_smooth = smoothdata(e(range2)/iter,'gaussian',10)
    ln_smooth = log(unif_EKK(range2))

    p=plot(log(unif_EKK(range2)),r_smooth,'-.'  )  ;
    p=plot(log(unif_EKK(range2)),e_smooth,'-.'  )  ;

    rho_interp = @(n) interp1(ln_smooth, r_smooth, log(n), 'pchip')  
    eps_interp = @(n) interp1(ln_smooth, e_smooth, log(n), 'pchip')  


%% From the Literature (EKK)
    rho     = rho_interp
    epsilon = eps_interp
    rho_bar = @(n)  (1/n).* integral(@(x)rho(x) ,0,n);
    ratio   = @(n) rho(n)./epsilon(n);
    C          = size(n_ij,1);
    warning('off', 'MATLAB:integral:MaxIntervalCountReached')
    top1 = arrayfun(@(lb,ub) integral(ratio,lb,ub),zeros(C),n_ij);
    bot1 = arrayfun(@(lb,ub) integral(rho,lb,ub)  ,zeros(C),n_ij)./arrayfun(epsilon,n_ij);
    precalculated_integrals = top1./bot1;
    FGFT = @(n_ii_hat)    GFT(n_ii_hat,X_ij,epsilon,rho,rho_bar,ratio,n_ij,x_ij,1,precalculated_integrals,'array',1,0);

    [n_ii_hat_EKK,~,~,~,~] = fsolve(FGFT,ones(size(n_ij,1),1)+1,optimoptions('fsolve','Display','iter'));
    [~,N_hat_EKK] = FGFT(n_ii_hat_EKK);
    elasticity       = -1/(sigma-1);
    W_hat_splineEKK_1 = diag(x_ij).^elasticity;
    W_hat_splineEKK_2 = (n_ii_hat_EKK.*N_hat_EKK).^elasticity;
    W_hat_splineEKK_3 = ( arrayfun(rho_bar,n_ii_hat_EKK.*diag(n_ij))./arrayfun(rho_bar,diag(n_ij))    ).^elasticity;
    W_hat_splineEKK   = W_hat_splineEKK_1 .* W_hat_splineEKK_2 .* W_hat_splineEKK_3;
    EKK_check = mean(abs(FGFT(n_ii_hat_EKK)))
    
    results_EKK = [log(W_hat_splineEKK) log(W_hat_splineEKK_1) log(W_hat_splineEKK_2) log(W_hat_splineEKK_3)];

    %% Plot
    w_term_EKK = @(flag,i) results_EKK(flag,(i+1)) ;
    sumEKK_13 = @(type) w_term_EKK(type,1);
    sumEKK_45 = @(type) w_term_EKK(type,2) + w_term_EKK(type,3) ;
    metric = @(c) (w_term_1(c,0)-w_term_EKK(c,0))./w_term_EKK(c,0);
    %% Wrapper to Plot Aggregate Data
        all = poor | rich;
        a = scatter(ln_n_ij(rich),metric(rich),[],[0.6350 0.0780 0.1840],'filled','o');
        hold on
        b = scatter(ln_n_ij(poor),metric(poor),[],[0 0.4470 0.7410],'filled','square');
        text(ln_n_ij(all),metric(all),labels(all),'VerticalAlignment','top','HorizontalAlignment','right')
        hold off
        xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
        % xlim([-.05 .010]);ylim([-.05 .010]);
        f =gca;
        legend([a b ],{'Developed', 'Developing'},'FontSize',14,'Location','best')
        f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};

    ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{EKK}})/{ \log \hat{W}_{EKK}}$'},'FontSize',16,'Interpreter','Latex');        
    exportgraphics(f,"../../Output/relative_12345_GFT_EKK.pdf")
    mean(abs(metric(all)))
% 
%     metric = @(c) (w_term_2(c,0)-w_term_EKK(c,0))./w_term_EKK(c,0);
%     GFT_agg_plot
%     ylabel({'$({\log \hat{W}_{Constant} - \log \hat{W}_{EKK}})/{ \log \hat{W}_{EKK}}$'},'FontSize',16,'Interpreter','Latex');        
%     ylim([-1 1]);
%     exportgraphics(f,"../../Output/crelative_12345_GFT_EKK.pdf")
% 
%     metric = @(c) (w_term_1(c,0)-w_term_EKK(c,0))./w_term_1(c,0);
%     GFT_agg_plot
%     ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{EKK}})/{ \log \hat{W}_{Semi}}$'},'FontSize',16,'Interpreter','Latex');        
%     ylim([-1 1]);
%     exportgraphics(f,"../../Output/crelative_12345_GFT_EKK_1.pdf")
%     
%     metric = @(c) (w_term_1(c,0)-w_term_EKK(c,0))./w_term_2(c,0);
%     GFT_agg_plot
%     ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{EKK}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');        
%     ylim([-1.2 1]);
%     exportgraphics(f,"../../Output/crelative_12345_GFT_EKK_2.pdf")
    