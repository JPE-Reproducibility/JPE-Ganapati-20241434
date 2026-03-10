function  Simulation_Data_v10

    data_dir = '../../Data/Int/WIOD_sampleB/';
    addpath('../GMM_estimation_v3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of Various Firm Models and Theta Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta           = 4;    %Simonovska  and  Waugh  (2014a)  estimate  a  trade  elasticity  of  4.10  or  4.27  depending  on  the  data  used.  
    sigma_lit       = 4;    %Bernard  et  al.  (2003).
    tau             = 1.83; % https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351 page 1128
                            % Melitz, M. J. and Redding, S. J. (2015). New trade models, new welfare implications. American Economic Review, 105(3):1105--46.
    B               = 1;


    % xx = -9.3:.01:-.1;
    % unif = exp(xx);

    xx = -9.3:.01:-.1;
    unif = exp(xx);


    % space           = .0001;
    % unif            = 0:space:(1);
    % unif_EKK        = 0:space:1;
    unif_EKK = unif;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Pareto (Data)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     z              = unif.^(-1/theta);
     e              = B.*(tau./z).^(1-sigma_lit);
     r              = e;

     x2 = 1:(size(unif,2)-1);

     % e_f             = @(n) B.*(tau./n).^(1-sigma)
     % r_int = arrayfun(@(n) integral(e_f,0,n),z)./z;
     % elast_r2_logn   = r(x2)./r_int(x2) -1


     elast_e_pareto = (log(e(1:(size(unif,2)-1)))-log(e(2:size(unif,2))))./(log(unif(1:(size(unif,2)-1)))-log(unif(2:size(unif,2))));
     elast_r_pareto = (log(r(1:(size(unif,2)-1)))-log(r(2:size(unif,2))))./(log(unif(1:(size(unif,2)-1)))-log(unif(2:size(unif,2))));

     rho_bar=arrayfun(@(l,u) integral(@(x) B.*(tau./(x.^(-1/theta))).^(1-sigma_lit) ,l,u),zeros(size(unif)),unif)./unif;
     elast_rbar_pareto = (log(rho_bar(1:(size(unif,2)-1)))-log(rho_bar(2:size(unif,2))))./(log(unif(1:(size(unif,2)-1)))-log(unif(2:size(unif,2))));


%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Log-Normal Data
% Fernandes, A., Klenow, P., Denisse Periola, M., Meleshchuk, S., and Rodriguez-Clare, A. (2017). The intensive margin in trade: Moving beyond pareto. Unpublished manuscript.
% https://siepr.stanford.edu/sites/default/files/publications/555wp_1.pdf


% find the bas et al / Head & Meyer & Thoneing - trade without pareto
% https://ideas.repec.org/a/eee/inecon/v108y2017icp1-19.html


% That paper doesn't estimate it, it simply uses Head Meyer AER where 0.797
% = sigma 
% We are going to just use that
%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_sigma         = .79;  
    n_mu            = 0;
    
    % Updated Code; Should work - need to check
    % Should this be 1-unif or just unif?
    z               = exp(n_mu + sqrt(2)*n_sigma*erfinv(2*(1-unif)-1));
    e               = B.*(tau./z).^(1-sigma_lit);
    %
    r               = e;

    e_f             = @(n) B.*(tau./n).^(1-sigma_lit);
    z_f             = @(n) exp(n_mu + sqrt(2)*n_sigma*erfinv(2*(1-n)-1));
    r_int = arrayfun(@(n) integral(@(n2) e_f(z_f(n2)),0,n),unif)./unif;
    % rr = @(n2) e_f(z_f(n2))

    % r_int_f = @(n) integral(@(n2) e_f(z_f(n2)),0,n)

    x2 = 1:(size(unif,2)-1);
    x1 = 2:size(unif,2);


%     elast_r2_logn    = (log(r(x2))-log(r(x1)))./(log(unif(x2))-log(unif(x1)));
    elast_e_logn    = (log(e(x2))-log(e(x1)))./(log(unif(x2))-log(unif(x1)));
 
    elast_r_logn   = r(x2)./r_int(x2) -1;
    
%     range = 1:1:900;
%     ln_ij = log(unif(range));
%     clf;
%          
%     p=plot(ln_ij,elast_e_logn(range),'-',...
%            ln_ij,elast_r_logn(range),'--');


    % rho_bar=arrayfun(@(l,u) integral(@(x) B.*(tau./(exp(n_mu + sqrt(2)*n_sigma*erfinv(2*x-1)))).^(1-sigma_lit) ,l,u),zeros(size(unif)),fliplr(unif))./unif;
    % elast_rbar_logn = (log(rho_bar(1:(size(unif,2)-1)))-log(rho_bar(2:size(unif,2))))./(log(unif(1:(size(unif,2)-1)))-log(unif(2:size(unif,2))));

    % scatter(r(x2(1:800)),r_int(1:800))


    %%%% FUTURE TO DO: CLEAN UP NOMENCALRE FOR r and r_int!!!!!!

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Truncated Pareto   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     L                  = 1;
     H                  = 2.85; % From https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351
                                % Melitz, M. J. and Redding, S. J. (2015). New trade models, new welfare implications. American Economic Review, 105(3):1105--46.
     theta_tp           = 4;    % Simonovska  and  Waugh  (2014a)  estimate  a  trade  elasticity  of  4.10  or  4.27  depending  on  the  data  used.  
     z                  = 1-(1-L.^(1/theta_tp).*unif.^(-1/theta_tp))./(1-(L/H).^(1/theta_tp));
     e                  = B.*(tau./z).^(1-sigma_lit);
     r                  = e;
     % elast_r_tpareto    = (log(r(1:(size(unif,2)-1)))-log(r(2:size(unif,2))))./(log(unif(1:(size(unif,2)-1)))-log(unif(2:size(unif,2))));
     elast_e_tpareto    = (log(e(1:(size(unif,2)-1)))-log(e(2:size(unif,2))))./(log(unif(1:(size(unif,2)-1)))-log(unif(2:size(unif,2))));
     


    e_f             = @(n) B.*(tau./n).^(1-sigma_lit);
    z_f             = @(n) 1-(1-L.^(1/theta_tp).*n.^(-1/theta_tp))./(1-(L/H).^(1/theta_tp));
    r_int = arrayfun(@(n) integral(@(n) e_f(z_f(n)),0,n),unif)./unif;     
    elast_r_tpareto   = r(x2)./r_int(x2) -1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Pareto with log-normal entry (EKK-lite)   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this only works with a uniform distribution (symmetric)

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
%     p=plot(ln_ij,smoothdata(elast_r_EKK(range),'gaussian',200),'-',...
%            ln_ij,smoothdata(elast_e_EKK(range),'gaussian',200),'--');


%% Plot Simulation results
    theta_pareto  = (1-sigma_lit).*(1+elast_r_pareto  - elast_e_pareto )./elast_e_pareto  ;
    theta_tpareto = (1-sigma_lit).*(1+elast_r_tpareto - elast_e_tpareto)./elast_e_tpareto ;
    theta_logn    = (1-sigma_lit).*(1+elast_r_logn    - elast_e_logn   )./elast_e_logn    ;
    theta_EKK     = (1-sigma_lit).*(1+elast_r_EKK(1:(size(unif_EKK,2)-1))     - elast_e_EKK(1:(size(unif_EKK,2)-1))    )./elast_e_EKK(1:(size(unif_EKK,2)-1))     ;

    %% For text
    FONTSIZE = 22;
    range = 1:1:900;
    range2 = 2:min(1400,size(unif_EKK,2)-1);
    ln_ij = log(unif(range));

    clf;
    p=plot(ln_ij,elast_e_pareto(range),'-',...
         ln_ij,elast_e_tpareto(range),'--',...
         ln_ij,elast_e_logn(range),':');
         % log(unif_EKK(range2)),smoothdata(elast_e_EKK(range2),'gaussian',10),'-.' ...
         % )  ;
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;  %  p(4).LineWidth = 2;  
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    ax.FontSize = 16; 
    legend({'Pareto','Truncated Pareto','Log Normal','EKK'},'FontSize',FONTSIZE)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('\partial ln \epsilon (n ) / \partial ln (n)','FontSize',FONTSIZE)
    xlim([-8 -1])
    legend('Location','best')
    saveas(gcf,  '../../Output/Epsilon_v8','epsc')
    
 
    clf;
    p=plot(ln_ij,elast_r_pareto(range),'-',...
         ln_ij,elast_r_tpareto(range),'--',...
         ln_ij,elast_r_logn(range),':');
         % log(unif_EKK(range2)),smoothdata(elast_r_EKK(range2),'gaussian',10),'-.');
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    ax.FontSize = 16; 
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;   % p(4).LineWidth = 2; 
    legend({'Pareto','Truncated Pareto','Log Normal','EKK'},'FontSize',FONTSIZE)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('\partial ln \rho (n ) / \partial ln (n)','FontSize',FONTSIZE)
    xlim([-8 -1])
    legend('Location','best')
    saveas(gcf,  '../../Output/Rho_v8','epsc')
    

    clf;
    p=plot(ln_ij,theta_pareto(range),'-',...
         ln_ij,theta_tpareto(range),'--',...
         ln_ij,theta_logn(range),':');
         % log(unif_EKK(range2)),smoothdata(theta_EKK(range2),'gaussian',2),'-.');
         p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;   % p(4).LineWidth = 2;  
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    ax.FontSize = 16; 
    legend({'Pareto','Truncated Pareto','Log Normal','EKK'},'FontSize',FONTSIZE)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
        ylabel('$\theta$','FontSize',24,'Interpreter','Latex')
    xlim([-8 -1])
    saveas(gcf,  '../../Output/Theta_v8','epsc')
    



    clf;
    p=plot(ln_ij,-1./elast_e_pareto(range),'-',...
         ln_ij,-1./elast_e_tpareto(range),'--',...
         ln_ij,-1./elast_e_logn(range));%,':',...
         % log(unif_EKK(range2)),-1./smoothdata(elast_e_EKK(range2),'gaussian',10),'-.')  ;
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;  %  p(4).LineWidth = 2;  
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    ax.FontSize = 16; 
    legend({'Pareto','Truncated Pareto','Log Normal','EKK'},'FontSize',FONTSIZE)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
        ylabel('$\theta^{e}$','FontSize',24,'Interpreter','Latex')
    xlim([-8 -1])
    legend('Location','best')
    saveas(gcf,  '../../Output/Extensive_v8','epsc')
    

    clf;
    p=plot(ln_ij,1-elast_r_pareto(range)./elast_e_pareto(range),'-',...
         ln_ij,1-elast_r_tpareto(range)./elast_e_tpareto(range),'--',...
         ln_ij,1-elast_r_logn(range)./elast_e_logn(range),':');
         % log(unif_EKK(range2)),1-smoothdata(elast_r_EKK(range2),'gaussian',10)./smoothdata(elast_e_EKK(range2),'gaussian',10),'-.')  ;
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;   % p(4).LineWidth = 2;  
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    ax.FontSize = 16; 
    legend({'Pareto','Truncated Pareto','Log Normal','EKK'},'FontSize',FONTSIZE)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('Intensive Margin (\theta^i)','FontSize',FONTSIZE)
    xlim([-8 -1])
    legend('Location','best')
    saveas(gcf,  '../../Output/Intensive_v8','epsc')
    

    clf;
    p=plot(ln_ij,elast_r_pareto(range),'-',...
         ln_ij,elast_r_tpareto(range),'--',...
         ln_ij,elast_r_logn(range),':');%,...
         % log(unif_EKK(range2)),smoothdata(elast_r_EKK(range2),'gaussian',10),'-.'    );
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;  %  p(4).LineWidth = 2;  
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    ax.FontSize = 16; 
    legend({'Pareto','Truncated Pareto','Log Normal','EKK'},'FontSize',FONTSIZE)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('$\theta^{c}$','FontSize',24,'Interpreter','Latex')
    xlim([-8 -1])
    legend('Location','best')
    saveas(gcf,  '../../Output/Thetac_v8','epsc')



    range1 = range;
    save 'literature' elast* ln_ij range1 unif_EKK range2 theta*


    return

%% Plot Simulation results

    kappa_tau = 1;
    kappa_f   = 0;
    sigma_lit     = 3.2;    %     sigma     = 3.9; % sigma = 2.66

    kappa_epsilon = 1/((sigma_lit-1)*kappa_tau+kappa_f);
    gamma_guess = [.2 0 0 0 0 .2 0 0 0 0];

    % Use Unweighted
    basefile = '2012_Tetiuw.csv'; name = '2012_T'; %    basefile = '2012_Tuw_gravity.csv'; name = '2012_T';


    %% ADD EXPLANATORY VARIABLES
    fix = 2;
    knots = 3;
    M_base = csvread(strcat(data_dir,basefile));
    gbase = make_data_combo(M_base,knots,'cubic','','',0,'all');    gbase.sigma = sigma_lit;    %   gbase.x = x;
    [~,gobase] =  GMM_wrapper_gravity( gbase, kappa_tau, kappa_epsilon,'base',fix,gamma_guess,'all' );
    start = size(gobase.ests2,2)-2*size(gobase.est_epsilon,2)-1;
    gamma_gravity = gobase.ests2((end-start):end);


%% Baseline (taking the first stage results as given)
    d1 = make_data_combo(M_base,1,'cubic','','',0,'all');
    d1.sigma = sigma_lit;  % d1.x = x;
    [~,o1] =  GMM_wrapper_gravity( d1, kappa_tau, kappa_epsilon,'base',0 ,gamma_gravity,'all');
    [ epsilon_elast_loglin,rho_elast_loglin,theta_loglin ] = make_elasticity_combo(ln_ij,d1.k,d1,o1,sigma_lit,'none');


    d4 = make_data_combo(M_base,3,'cubic','','',0,'all');    d4.sigma = sigma_lit;    %   d4.x = x;
    [~,o4] =  GMM_wrapper_gravity( d4, kappa_tau, kappa_epsilon,'base',0 ,gamma_gravity,'all');
    [ epsilon_elast_spline,rho_elast_spline,theta_spline ] = make_elasticity_combo(ln_ij,d4.k,d4,o4,sigma_lit,'none');


    %% Wealth Origin
    fix = 0;
    knots = 34;
    d2q_Wo = make_data_combo(M_base,knots,'cubic','wealth_oI','',0,'all'); 
    d2q_Wo.sigma = sigma_lit; % d2q_Wo.x = x;
    [~,o2q_Wo] =  GMM_wrapper_gravity( d2q_Wo, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');

    [ epsilon_elast_o0spline,rho_elast_o0spline,theta_o0spline ] = make_elasticity_combo(ln_ij,d2q_Wo.k,d2q_Wo,o2q_Wo,sigma_lit,'wealth_oI0');
    [ epsilon_elast_o1spline,rho_elast_o1spline,theta_o1spline ] = make_elasticity_combo(ln_ij,d2q_Wo.k,d2q_Wo,o2q_Wo,sigma_lit,'wealth_oI1');

    clf;
    p=plot(ln_ij,theta_pareto(range),'-',...
        ln_ij,theta_tpareto(range),'--',...
        ln_ij,theta_logn(range),':',...
        ln_ij,theta_loglin,':',...
        ln_ij,theta_spline,'-.',...
        ln_ij,theta_o0spline,'-.' , ...
        ln_ij,theta_o1spline,'-.' );
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;    p(4).LineWidth = 3;    p(5).LineWidth = 2;  p(6).LineWidth = 2;  p(7).LineWidth = 2;
    ax              = gca; p(4).Color = 'black';
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    legend({'Pareto','Truncated Pareto','Log Normal','Constant Elasticity','Semiparametric','Poor Origin','Rich Origin'},'FontSize',16)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('\theta_{ij}','FontSize',FONTSIZE)
        set(gca, 'YScale', 'log')

    saveas(gcf,  '../../Output/Theta_overlay_v8','pdf')
    saveas(gcf,  '../../Output/Theta_overlay_v8','epsc')
   


    clf;
    p=plot(ln_ij,theta_pareto(range),'-',...
        ln_ij,theta_tpareto(range),'--',...
        ln_ij,theta_logn(range),':',...
        ln_ij,theta_loglin,':',...
        ln_ij,theta_spline,'-.');
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;    p(4).LineWidth = 4;    p(5).LineWidth = 4; 
    ax              = gca; p(5).Color = 'black'; p(5).LineStyle= '-'; p(4).Color = "#777777";p(5).LineStyle= '-'; 
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    legend({'Pareto','Truncated Pareto','Log Normal','Constant Elasticity','Semiparametric','Poor Origin','Rich Origin'},'FontSize',16)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('\theta_{ij}','FontSize',FONTSIZE)
    set(gca, 'YScale', 'log')
    saveas(gcf,  '../../Output/Theta_overlay2_v8','pdf')
    saveas(gcf,  '../../Output/Theta_overlay2_v8','epsc')
   

    p=plot(ln_ij,elast_r_pareto(range),'-',...
     ln_ij,elast_r_tpareto(range),'--',...
     ln_ij,elast_r_logn(range),':',...
     ln_ij,rho_elast_loglin,':',...
     ln_ij,rho_elast_spline,'-.',...
        ln_ij,rho_elast_o0spline,'-.' , ...
        ln_ij,rho_elast_o1spline,'-.' );
    % ylim([-2.5 0 ]);
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;    p(4).LineWidth = 2;    p(5).LineWidth = 2;  p(6).LineWidth = 2;  p(7).LineWidth = 2;
    legend({'Pareto','Truncated Pareto','Log Normal','Constant Elasticity','Semiparametric','Poor Origin','Rich Origin'},'FontSize',16)
    legend('boxoff')
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('\partial ln \rho (n ) / \partial ln (n)','FontSize',FONTSIZE)
    legend('Location','best')
    saveas(gcf,  '../../Output/Rho_overlay_v8','epsc')
    saveas(gcf,  '../../Output/Rho_overlay_v8','pdf')

    

    p=plot(ln_ij,elast_e_pareto(range),'-',...
         ln_ij,elast_e_tpareto(range),'--',...
         ln_ij,elast_e_logn(range),':',...
         ln_ij,epsilon_elast_loglin,':',...
         ln_ij,epsilon_elast_spline,'-.',...
        ln_ij,epsilon_elast_o0spline,'-.' , ...
        ln_ij,epsilon_elast_o1spline,'-.' );
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;    p(4).LineWidth = 2;    p(5).LineWidth = 2;  p(6).LineWidth = 2;  p(7).LineWidth = 2;
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    legend({'Pareto','Truncated Pareto','Log Normal','Constant Elasticity','Semiparametric','Poor Origin','Rich Origin'},'FontSize',16)
    legend('boxoff')
    ylim([-3 0 ]);
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('\partial ln \epsilon (n ) / \partial ln (n)','FontSize',FONTSIZE)
    legend('Location','best')
    saveas(gcf,  '../../Output/Epsilon_overlay_v8','epsc')
    saveas(gcf,  '../../Output/Epsilon_overlay_v8','pdf')


    p=plot(ln_ij,-1./elast_e_pareto(range),'-',...
         ln_ij,-1./elast_e_tpareto(range),'--',...
         ln_ij,-1./elast_e_logn(range),':',...
         ln_ij,-1./epsilon_elast_loglin,':',...
         ln_ij,-1./epsilon_elast_spline,'-.',...
        ln_ij,-1./epsilon_elast_o0spline,'-.' , ...
        ln_ij,-1./epsilon_elast_o1spline,'-.' );
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;    p(4).LineWidth = 2;    p(5).LineWidth = 2;  p(6).LineWidth = 2;  p(7).LineWidth = 2;
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    legend({'Pareto','Truncated Pareto','Log Normal','Constant Elasticity','Semiparametric','Poor Origin','Rich Origin'},'FontSize',16)
    legend('boxoff')
%     ylim([-3 0 ]);
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('Extensive Margin','FontSize',FONTSIZE)
    legend('Location','best')
    saveas(gcf,'../../Output/Extensive_overlay_v8','epsc')
    saveas(gcf, '../../Output/Extensive_overlay_v8','pdf')
    
    






    p=plot(ln_ij,1-elast_r_pareto(range)./elast_e_pareto(range),'-',...
         ln_ij,1-elast_r_tpareto(range)./elast_e_tpareto(range),'--',...
         ln_ij,1-elast_r_logn(range)./elast_e_logn(range),':',...
         ln_ij,1-rho_elast_loglin./epsilon_elast_loglin,':',...
         ln_ij,1-rho_elast_spline./epsilon_elast_spline,'-.',...
        ln_ij,1-rho_elast_o0spline./epsilon_elast_o0spline,'-.' , ...
        ln_ij,1-rho_elast_o1spline./epsilon_elast_o1spline,'-.' );
    p(1).LineWidth = 2;    p(2).LineWidth = 2;    p(3).LineWidth = 2;    p(4).LineWidth = 2;    p(5).LineWidth = 2;  p(6).LineWidth = 2;  p(7).LineWidth = 2;
    ax              = gca;
    ax.XTick   = [-9.210340372 -6.907755279 -5.298317367 -4.605170186 -2.995732274 -2.302585093 -0.693147181 0 ];
    ax.XTickLabel = {'0.01%' '0.1%','0.5%','1%','5%','10%','50%','100%'};
    legend({'Pareto','Truncated Pareto','Log Normal','Constant Elasticity','Semiparametric','Poor Origin','Rich Origin'},'FontSize',16)
    legend('boxoff')
%     ylim([-3 0 ]);
    xlabel('Log Exporter Firm Share','FontSize',FONTSIZE)
    ylabel('Intensive Margin','FontSize',FONTSIZE)
    legend('Location','best')
    saveas(gcf,'../../Output/Intensive_overlay_v8','epsc')
    saveas(gcf,'../../Output/Intensive_overlay_v8','pdf')



end