 % Data Needs
    clear

%% Using Estimates
    year = '2012B';    
    addpath( '../GMM_estimation_v3')
    data_dir  =  '../../Data/Int/WIOD_sampleB/';
    excel_output = [ '../../Output/GFT_results_revision_' year '.xlsx'];
    basefile = '2012_Tetiuw.csv';


    n_ij            = csvread([ data_dir 'n_ij_'      year '.csv']);
    x_ij            = csvread([ data_dir 'x_ij_'      year '.csv']);
    X_ij            = csvread([ data_dir 'XX_ij_'     year '.csv']);
    Agg_i           = csvread([ data_dir 'balance_i_' year '.csv']);
    kappa_i         = Agg_i(:,3);
    labels          = textread([data_dir 'l_i_'       year '.csv'],'%q');
    G_ij            = csvread([ data_dir 'G_ij_'      year '.csv']);
    G4_ij           = csvread([ data_dir 'G4_ij_'     year '.csv']);
    sample_selection= csvread([ data_dir 'Sample_'    year '.csv']);

    kappa_tau = 1;
    kappa_f   = 0;
    sigma     = 3.2; 
    kappa_epsilon = 1/((sigma-1)*kappa_tau+kappa_f);
    % gamma_gravity = [ 0.3920   -0.0408   -0.0131   -0.2596   -0.1635    0.3402   -0.1012    0.0039   -0.0985   -0.0671];
    % gamma_gravity =    [0.3520   -0.0214   -0.0189   -0.2324   -0.1502    0.2812   -0.0834   -0.0015   -0.0576   -0.0471];
    load( '../GMM_estimation_v3/gamma_gravity.mat');

    M_base = csvread(strcat(data_dir,basefile));
    
    C               = size(n_ij,1);
    y_ij = X_ij ./ repmat(sum(X_ij,2),1,C);
    x_ij = X_ij ./ repmat(sum(X_ij,1),C,1);


%% Using Estimates - Constant Elasticity
    fix = 0; 
    knots = 1;
    d1 = make_data_combo(M_base,knots,'cubic','','',0,'all');     d1.sigma = sigma;    
    [~,o1] =  GMM_wrapper_gravity( d1, kappa_tau, kappa_epsilon,'base',0,gamma_gravity,'all' );
    rho_power       = o1.est_rho;  
    e_power         = o1.est_epsilon;  

    %convert to macros
    epsilon_bar     = @(n) n.^e_power;
    rho             = @(n) (rho_power+1)*n.^rho_power;
    rho_bar         = @(n) n.^(rho_power);
    ratio           = @(n) rho(n)./epsilon_bar(n);

    % we need to start epsilon away from zero
    warning('on','MATLAB:integral:MinStepSize')

    top1 = arrayfun(@(lb,ub) integral(ratio,lb,ub),zeros(C),n_ij);
    bot1 = arrayfun(@(lb,ub) integral(rho,lb,ub)  ,zeros(C),n_ij)./epsilon_bar(n_ij);
    precalculated_integrals = top1./bot1; % Also Equal to (rho_power+1)/(rho_power+1-e_power);

    % Rho or Rho_bar?
    FGFT = @(n_ii_hat) GFT_v9(n_ii_hat,X_ij,epsilon_bar,rho_bar,rho_bar,ratio,n_ij,x_ij,1,precalculated_integrals,'noarray',0,0);
    [n_ii_hat_generalized_pareto,~,~,~,~] = fsolve(FGFT,ones(size(n_ij,1),1)+1,optimoptions('fsolve','Display','iter'));
    [~,N_hat_generalized_pareto] = FGFT(n_ii_hat_generalized_pareto);

    elasticity = -1/(sigma-1);
    W_hat_gen_pareto_1 = diag(x_ij).^elasticity;
    W_hat_gen_pareto_2 = (N_hat_generalized_pareto).^elasticity;
    W_hat_gen_pareto_3 = (n_ii_hat_generalized_pareto.*rho_bar(n_ii_hat_generalized_pareto.*diag(n_ij))./rho_bar(diag(n_ij))).^elasticity;
    W_hat_gen_pareto   =  W_hat_gen_pareto_1 .*  W_hat_gen_pareto_2 .* W_hat_gen_pareto_3;


    share = mean(n_ij-diag(diag(n_ij)),2);
    Results_main = [  diag(x_ij) share log(W_hat_gen_pareto_1) log(W_hat_gen_pareto_2)  log(W_hat_gen_pareto_3) log(W_hat_gen_pareto) diag(G_ij)/3  N_hat_generalized_pareto n_ii_hat_generalized_pareto];
    T = array2table(Results_main,'VariableNames', {'x_ii' ,'Mean Exporter Share', 'Term1' 'Term2','Term3','Welfare','Developing Status','N_hat','n_ii_hat'});
    T = [labels T];
    writetable(T,excel_output,'sheet','Generalized Pareto');


    % Alternative Calculation
    W_hat_gen_pareto_alt = epsilon_bar((n_ii_hat_generalized_pareto.*N_hat_generalized_pareto)).^elasticity;

    % Check to make sure the math works
    elasticity_pareto = -e_power/(1+rho_power-e_power)/(1-sigma);
    W_hat_pareto_simple = diag(x_ij).^elasticity_pareto;


    % This should all be zero!
    fprintf( 'ERROR: % 10.5f\t\n',    mean(abs(W_hat_gen_pareto./W_hat_pareto_simple))-1)
    fprintf( 'ERROR: % 10.5f\t\n',    mean(abs(W_hat_pareto_simple./W_hat_gen_pareto_alt))-1)

    % This is a measure of the rounding error!!
    for i = 1:max(size(W_hat_pareto_simple))
        fprintf('%s\t% 10.5f\t% 10.5f\t% 10.5f\t% 10.5f\n',labels{i},100*log(W_hat_pareto_simple(i)),100*log(W_hat_gen_pareto_alt(i)),100*log(W_hat_gen_pareto(i)),100*(log(W_hat_gen_pareto(i))-log(W_hat_pareto_simple(i)))/log(W_hat_pareto_simple(i))) ;
    end

    

%% %%%%%%%%%%%%
% Groups - origin_destination groups
    knots = 3;
    dGA = make_data_combo(M_base,knots,'cubic','wealth_odI','',0,'all'); 
    dGA.sigma = sigma;
    [~,oGA] =  GMM_wrapper_gravity( dGA, kappa_tau, kappa_epsilon,'base',0,gamma_gravity,'all' );

    GA_ij = G_ij;

    % Get Functions
    GProjection         = @(n,G) evknots_makeG(dGA.kG,log(n),G,2) ;
    GDerivitive         = @(n,G) evknots_make_derivG(dGA.kG,log(n),G,2) ;
    GrhoA               = @(n,G) exp((GProjection(n,G)*oGA.est_rho')') ;
    GepsilonA           = @(n,G) exp((GProjection(n,G)*oGA.est_epsilon')');
    Gelast_epsilon      = @(x,G) reshape( GDerivitive((x(:)),G(:))*oGA.est_epsilon',1,[]);
    Gelast_rho          = @(x,G) reshape( GDerivitive((x(:)),G(:))*oGA.est_rho',1,[]);
    Gtheta              = @(n,G) (1-sigma)*(1+Gelast_rho(n,G) - Gelast_epsilon(n,G))./Gelast_epsilon(n,G);
    GratioA             = @(n,G) GrhoA(n,G)./GepsilonA(n,G).*(1+Gelast_rho(n,G));
    Grho_barA           = @(n,G) GrhoA(n,G).*(1+Gelast_rho(n,G));

    G_values = unique(unique(GA_ij));
    gamma_ij = zeros(C);

    for g = 1:size(G_values,1)
        gg = G_values(g);
        GratioA0 = @(nn) GratioA(nn,gg);              
        GrhoA0 = @(nn) GrhoA(nn,gg);
        GepsilonA0 = @(nn) GepsilonA(nn,gg); 
        gamma0_ij = arrayfun(@(n) GepsilonA0(n)./n./GrhoA0(n),n_ij).*arrayfun(@(lb,ub) integral(GratioA0,lb,ub) , zeros(C),n_ij);
        gamma_ij(GA_ij==gg) = gamma0_ij(GA_ij==gg);
    end

    % Equation eq:Ni_GT
    denominator = 1-sum(y_ij.*gamma_ij,2); 

    % Don't Fix N_i
    FGFT = @(n_ii_hat)  GFTG_v10(n_ii_hat,n_ij,x_ij,GA_ij,denominator,GrhoA,GepsilonA,GratioA,Grho_barA,1);
    [n_ii_hat_od,~,~,~,~] = fsolve(FGFT,ones(size(n_ij,1),1),optimoptions('fsolve','Display','iter','OptimalityTolerance',1e-10));
    [~,N_hat_od] = FGFT(n_ii_hat_od);
    W_hat_spline_od = (diag(x_ij).*n_ii_hat_od.*N_hat_od.*( arrayfun(GrhoA,n_ii_hat_od.*diag(n_ij),diag(GA_ij))./arrayfun(GrhoA,diag(n_ij),diag(GA_ij))    )).^(1/(1-sigma));

    % Average Exporter Firm Share
    share = mean(n_ij-diag(diag(n_ij)),2);

    % Save Results
    elasticity      = 1/(1-sigma);
    W_hat_spline_od_1 = diag(x_ij).^elasticity;
    W_hat_spline_od_2 = (N_hat_od).^elasticity;
    W_hat_spline_od_3 = ( n_ii_hat_od.*( arrayfun(GrhoA,n_ii_hat_od.*diag(n_ij),diag(GA_ij))./arrayfun(GrhoA,diag(n_ij),diag(GA_ij))    )   ).^elasticity;
    W_hat_spline_od_c = W_hat_spline_od_1 .* W_hat_spline_od_2 .* W_hat_spline_od_3;
    Results_main = [  diag(x_ij) share log(W_hat_spline_od_1) log(W_hat_spline_od_2)  log(W_hat_spline_od_3) log(W_hat_spline_od_c) diag(G_ij)/3 N_hat_od n_ii_hat_od];
    T = array2table(Results_main,'VariableNames', {'x_ii' ,'Mean Exporter Share', 'Term1' 'Term2','Term3','Welfare','Developing Status','N_hat','n_ii_hat'});
    T = [labels T];
    writetable(T,excel_output,'sheet','SplineOD');

%% Decomposition
    weights     = sum(X_ij,1)';
    weights     = weights/sum(weights);
    
    results_linear = [log(W_hat_gen_pareto) log(W_hat_gen_pareto_1) log(W_hat_gen_pareto_2) log(W_hat_gen_pareto_3)];
    results_spl_od = [log(W_hat_spline_od) log(W_hat_spline_od_1) log(W_hat_spline_od_2) log(W_hat_spline_od_3)];

    OD      = diag(GA_ij);
    x_ii    = diag(x_ij);
    rich    = (OD==0) & (sample_selection == 3);
    poor    = (OD==3) & (sample_selection == 3);
    share   = mean(n_ij-diag(diag(n_ij)),2);
    ln_n_ij = log(share);
    share = mean(n_ij-diag(diag(n_ij)),2);
    ln_n_ij = log(share);

    % Decomposition Exact - Weighted
    fprintf('\n\nDecomposition Exact - Weighted\n')
    fprintf('Constant Elasticity \t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',100*sum(results_linear.*weights))
    fprintf('  Developed Origins \t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',100*sum(rich.*results_linear.*weights)/sum(rich.*weights))
    fprintf('  Developing Origins\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',100*sum(poor.*results_linear.*weights)/sum(poor.*weights))
    fprintf('Spline 4-Group      \t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',100*sum(results_spl_od.*weights))
    fprintf('  Developed Origins \t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',100*sum(rich.*results_spl_od.*weights)/sum(rich.*weights))
    fprintf('  Developing Origins\t %10.3f\t%10.3f\t%10.3f\t%10.3f\n',100*sum(poor.*results_spl_od.*weights)/sum(poor.*weights))

    w_term_1 = @(flag,i) results_spl_od(flag,(i+1)) ;
    w_term_2 = @(flag,i) results_linear(flag,(i+1)) ;
    sum1_13 = @(type) w_term_1(type,1);
    sum2_13 = @(type) w_term_2(type,1);
    sum1_45 = @(type) w_term_1(type,2) + w_term_1(type,3) ;
    sum2_45 = @(type) w_term_2(type,2) + w_term_2(type,3) ;

    %%%%%%%%%%%%%%%%%%%%%%
    %% Plot the Levels  %%
    %%%%%%%%%%%%%%%%%%%%%%
    %% Panel 2 N_i
        xmetric  = @(c) log(1./(x_ii(c)));
        ymetric1 = @(c) log(1./(N_hat_od(c)) );
        ymetric2 = @(c) log(1./(N_hat_generalized_pareto(c)) );
        GFT_division
        xlabel('Import Share $\log(1/x_{ii})$','FontSize',16,'Interpreter','Latex');            
        ylabel('$\log \hat{N}_{i}$','FontSize',16,'Interpreter','Latex')
        exportgraphics(f,"../../Output/2_groups_x_Ni.pdf")
    
        xmetric  = @(c) ln_n_ij(c);
        GFT_division
        f =gca;f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
        xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
        ylabel('$\log \hat{N}_{i}$','FontSize',16,'Interpreter','Latex')
        exportgraphics(f,"../../Output/2_groups_n_Ni.pdf")


    %% Panel 2b - 
    % Profit Shares

    gamma_new_ij = zeros(C);

    for g = 1:size(G_values,1)
        gg = G_values(g);
        GratioA0 = @(nn) GratioA(nn,gg);              
        GrhoA0 = @(nn) GrhoA(nn,gg);
        GepsilonA0 = @(nn) GepsilonA(nn,gg); 
        gamma0_ij = arrayfun(@(n) GepsilonA0(n)./n./GrhoA0(n),n_ij).*arrayfun(@(lb,ub) integral(GratioA0,lb,ub) , zeros(C),n_ij.*(n_ii_hat_od));
        gamma_new_ij(GA_ij==gg) = gamma0_ij(GA_ij==gg);
    end
    
    y_new_ij = eye(C);
    

    s_i_pi = 1-sum((y_ij.*gamma_ij),2);
    s_i_pi_new = 1-sum((y_new_ij.*gamma_new_ij),2);
    
    
    s_i_piD     = diag(y_ij).*(1-diag(gamma_ij)) ./ s_i_pi
    s_i_piD_new = diag(y_new_ij).*(1-diag(gamma_new_ij)) ./ s_i_pi_new
    
    hat_s_i_piD = -log(s_i_piD_new./s_i_piD)
    ymetric1 = @(c) ((hat_s_i_piD(c)) );
    ymetric2 = @(c) ((hat_s_i_piD(c)) );
    xmetric  = @(c) ln_n_ij(c);
    GFT_division2   
    legend([a b ],{'Semiparametric - Developed', 'Semiparametric - Developing' },'FontSize',14,'Location','southwest')
    f =gca;f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
%     ylabel('Change in Domestic Share of Profits: $\log \hat{s}_{i}^{\pi,D}$','FontSize',16,'Interpreter','Latex');
    ylabel('Change in Domestic Share of Profits','FontSize',16,'Interpreter','Latex');
    exportgraphics(f,"../../Output/2_groups_p_nii.pdf")
    
    
    ymetric1 = @(c) -log((N_hat_od(c)) );
    ymetric2 = @(c) -log((N_hat_od(c)) );
    GFT_division2   
    legend([a b ],{'Semiparametric - Developed', 'Semiparametric - Developing' },'FontSize',14,'Location','northeast')
    f =gca;f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
%     ylabel('Change in Profit Share: $\log \hat{s}_{i}^{\pi}$','FontSize',16,'Interpreter','Latex');
    ylabel('Change in Profit Share','FontSize',16,'Interpreter','Latex');
    exportgraphics(f,"../../Output/2_groups_p2_nii.pdf")
    
    
    
    ymetric1 = @(c) -log((N_hat_od(c)) ) + ((hat_s_i_piD(c)) );
    ymetric2 = @(c) -log((N_hat_od(c)) ) + ((hat_s_i_piD(c)) );
    GFT_division2   
    legend([a b ],{'Semiparametric - Developed', 'Semiparametric - Developing' },'FontSize',14,'Location','southwest')
    f =gca;f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
%     ylabel('Change in Domestic Profit Share of Income: $\log \hat{s}_{i}^{\pi}$','FontSize',16,'Interpreter','Latex');
    ylabel('Change in Domestic Profit Share of Income','FontSize',16,'Interpreter','Latex');
    exportgraphics(f,"../../Output/2_groups_p3_nii.pdf")
    
    
    histogram(diag(x_ij)./s_i_pi)
    xlabel('$x_{ii}/s_{i}^\pi$','FontSize',20,'Interpreter','Latex')

    histogram((sigma-1)/sigma*diag(x_ij)./s_i_pi)

    Y = (sigma-1)/sigma*diag(x_ij)./s_i_pi
    ymetric1 = @(c) -((Y(c)) );
    ymetric2 = @(c) -((Y(c)) );
    GFT_division2   
    legend([a b ],{'Semiparametric - Developed', 'Semiparametric - Developing' },'FontSize',14,'Location','northwest')
    f =gca;f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
    ylabel('$x_{ii}/s_{i}^\pi$','FontSize',20,'Interpreter','Latex')


    rho = (sigma-1)/sigma*diag(x_ij)./s_i_pi ./ s_i_piD
    histogram(rho.*rich)
    hold on
    histogram(rho.*poor)
    hold off
    xlabel('$\rho$','FontSize',20,'Interpreter','Latex')



    edges = linspace(.5, 3.5, 15);
    histogram(rho(rich), 'BinEdges', edges, 'FaceAlpha', 0.5, 'DisplayName', 'Developed Countries');
    hold on;
    histogram(rho(poor), 'BinEdges', edges, 'FaceAlpha', 0.5, 'DisplayName', 'Developing Countries');
    xlabel('$\rho$','FontSize',20,'Interpreter','Latex')
    legend('show', 'FontSize', 14,'Location','northwest');
    hold off;
    f = gcf;
    exportgraphics(f,"../../Output/GFT_rho.pdf")

    edges = linspace(.5, 3.5, 15);
    histogram(rho, 'BinEdges', edges, 'FaceAlpha', 0.5, 'DisplayName', 'Developed Countries');
    % hold on;
    % histogram(rho(poor), 'BinEdges', edges, 'FaceAlpha', 0.5, 'DisplayName', 'Developing Countries');
    xlabel('$\rho$','FontSize',20,'Interpreter','Latex')
    % legend('show', 'FontSize', 14);
    f = gcf;
    exportgraphics(f,"../../Output/GFT_rho_one.pdf")



    median(rho(rich))
    median(rho(poor))
    median(rho)

    
    % % --- Assume 'rho', 'rich', and 'poor' exist in your workspace ---
    % 
    % % 1. Calculate the median values
    % median_rich = median(rho(rich));
    % median_poor = median(rho(poor));
    % 
    % txt_rich = sprintf('Median = %.1f', median_rich);
    % txt_poor = sprintf('Median = %.1f', median_poor);
    % 
    % edges = linspace(0, 4, 20);
    % f = figure; % Create a new figure window AND capture its handle in 'f'
    % 
    % h1 = histogram(rho(rich), 'BinEdges', edges, 'FaceAlpha', 0.5, 'DisplayName', 'Developed Countries');
    % hold on;
    % 
    % h2 = histogram(rho(poor), 'BinEdges', edges, 'FaceAlpha', 0.5, 'DisplayName', 'Developing Countries');
    % 
    % xline(median_rich, 'LineWidth', 2, 'LineStyle', '--', ...
    %       'Color', 'black', 'HandleVisibility', 'off');
    %       
    % xline(median_poor, 'LineWidth', 2, 'LineStyle', ':', ...
    %       'Color', 'black', 'HandleVisibility', 'off');
    % 
    % yl = ylim; % Get the top of the y-axis
    % text(median_rich, yl(2) * 0.9, txt_rich, ...
    %      'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');
    %      
    % text(median_poor, yl(2) * 0.8, txt_poor, ...
    %      'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    %      'FontSize', 12, 'Color', 'black', 'FontWeight', 'bold');
    % 
    % hold off; % Release the plot hold
    %     
    %     xlabel('$\rho$','FontSize',20,'Interpreter','Latex')
    %     
    %     % This legend will now only show 'Developed Countries' and 'Developing Countries'
    %     legend('show', 'FontSize', 14, 'Location', 'best'); 
    %     
    %     exportgraphics(f,"../../Output/GFT_rho.pdf")






    rho = (sigma-1)/sigma*diag(x_ij)./s_i_pi ./ s_i_piD
    ymetric1 = @(c) ((rho(c)) );
    ymetric2 = @(c) ((rho(c)) );
    GFT_division2   
    legend([a b ],{'Semiparametric - Developed', 'Semiparametric - Developing' },'FontSize',14,'Location','northwest')
    f =gca;f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
    ylabel('$\rho$','FontSize',20,'Interpreter','Latex')

    Y = (sigma-1)/sigma*diag(x_ij)./s_i_pi.*log(W_hat_spline_od)
    ylabel('$x_{ii}/s_{i}^\pi\times d \log(\hat{W})$','FontSize',20,'Interpreter','Latex')



    %% Panel 3 - n_ii
        xmetric  = @(c) log(1./(x_ii(c)));
        ymetric1 = @(c) log(1./(n_ii_hat_od(c)) );
        ymetric2 = @(c) log(1./(n_ii_hat_generalized_pareto(c)) );
        GFT_division   
        xlabel('Import Share $\log(1/x_{ii})$','FontSize',16,'Interpreter','Latex');            
        ylabel('$\log \hat{n}_{ii}$','FontSize',16,'Interpreter','Latex');
        exportgraphics(f,"../../Output/2_groups_x_nii.pdf")

        xmetric  = @(c) ln_n_ij(c);
        GFT_division
        f =gca;f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
        xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
        ylabel('$\log \hat{n}_{ii}$','FontSize',16,'Interpreter','Latex');
        exportgraphics(f,"../../Output/2_groups_n_nii.pdf")

    %% Output Plots For 1-2-3-4-5
        metric = @(c) (w_term_1(c,0)-w_term_2(c,0))./w_term_2(c,0);
        GFT_agg_plot
        ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
        ylim([-.5 .6]);

        exportgraphics(f,"../../Output/relative_12345_GFT.pdf")

    %% Text Points
        rich_deviation = mean((metric(rich)))
        poor_deviation = mean((metric(poor)))
        min((metric(poor)))
        [i1,i2] = min((metric(poor)));ll = labels(poor);ll(i2)
        
         max((metric(rich)))
        [i1,i2] = max((metric(rich)));ll = labels(rich);ll(i2)
          


    %% Output Plots For 1-2-3
        metric = @(c) (sum1_13(c)-sum2_13(c))./w_term_2(c,0);
        GFT_agg_plot
        ylabel({'$({\log W^{Neoclassical}_{Semi} - \log W^{Neoclassical}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f, "../../Output/relative_123_GFT.pdf")

    %% Output Plots For4-5
        metric = @(c) (sum1_45(c)-sum2_45(c))./w_term_2(c,0);
        GFT_agg_plot 
        ylabel({'$({\log W^{Firm}_{Semi} - \log W^{Firm}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,"../../Output/relative_45_GFT.pdf")

    %% Output Plots For 4
        metric = @(c) (w_term_1(c,2)-w_term_2(c,2))./w_term_2(c,0);
        GFT_agg_plot 
        ylabel({'$({\log W^{Entry}_{Semi} - \log W^{Entry}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex'); 
        ylim([-.4 1.6]);
        exportgraphics(f,"../../Output/relative_4_GFT.pdf")

    %% Output Plots For 5
        metric = @(c) (w_term_1(c,3)-w_term_2(c,3))./w_term_2(c,0);
        GFT_agg_plot
        ylabel({'$({\log W^{Selection}_{Semi} - \log W^{Selection}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
        ylim([-1.4 .6]);
        exportgraphics(f, "../../Output/relative_5_GFT.pdf")


%% From the Literature (TP)
    tau        = 1.83; % https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351 page 1128
                       % Melitz, M. J. and Redding, S. J. (2015). New trade models, new welfare implications. American Economic Review, 105(3):1105--46.
    B          = 1;
    L          = 1;
    H          = 2.85; % From https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.20130351
                            % Melitz, M. J. and Redding, S. J. (2015). New trade models, new welfare implications. American Economic Review, 105(3):1105--46.
    theta_tp   = 4;    % Simonovska  and  Waugh  (2014a)  estimate  a  trade  elasticity  of  4.10  or  4.27  depending  on  the  data  used.
    z_TP       = @(n) 1-(1-L.^(1/theta_tp).*n.^(-1/theta_tp))./(1-(L/H).^(1/theta_tp));
    rho_TP     = @(n)  B.*(tau./z_TP(n)).^(1-sigma);
    epsilon_TP = @(n)  rho_TP(n);
    rho_bar_TP = @(n)  (1/n).* integral(@(x)B.*(tau./( 1-(1-L.^(1/theta_tp).*x.^(-1/theta_tp))./(1-(L/H).^(1/theta_tp)))).^(1-sigma) ,0,n);
    ratio_TP   = @(n) rho_TP(n)./epsilon_TP(n);
    C          = size(n_ij,1);

    top1 = arrayfun(@(lb,ub) integral(ratio_TP,lb,ub),zeros(C),n_ij);
    bot1 = arrayfun(@(lb,ub) integral(rho_TP,lb,ub)  ,zeros(C),n_ij)./arrayfun(epsilon_TP,n_ij);
    precalculated_integrals = top1./bot1;
    FGFT = @(n_ii_hat)    GFT(n_ii_hat,X_ij,epsilon_TP,rho_TP,rho_bar_TP,ratio_TP,n_ij,x_ij,1,precalculated_integrals,'array',1,0);
    % FGFT = @(n_ii_hat) GFT_v9(n_ii_hat,X_ij,epsilon_TP,rho_TP,rho_bar_TP,ratio_TP,n_ij,x_ij,1,precalculated_integrals,'noarray',0,0);

    [n_ii_hat_TP,~,~,~,~] = fsolve(FGFT,ones(size(n_ij,1),1)+1,optimoptions('fsolve','Display','iter'));
    [~,N_hat_TP] = FGFT(n_ii_hat_TP);
    elasticity       = -1/(sigma-1);
    W_hat_splineTP_1 = diag(x_ij).^elasticity;
    W_hat_splineTP_2 = (n_ii_hat_TP.*N_hat_TP).^elasticity;
    W_hat_splineTP_3 = ( arrayfun(rho_bar_TP,n_ii_hat_TP.*diag(n_ij))./arrayfun(rho_bar_TP,diag(n_ij))    ).^elasticity;
    W_hat_splineTP   = W_hat_splineTP_1 .* W_hat_splineTP_2 .* W_hat_splineTP_3;
    TP_check = mean(abs(FGFT(n_ii_hat_TP)))
    
    results_TP = [log(W_hat_splineTP) log(W_hat_splineTP_1) log(W_hat_splineTP_2) log(W_hat_splineTP_3)];

    %% Plot
    w_term_TP = @(flag,i) results_TP(flag,(i+1)) ;
    sumTP_13 = @(type) w_term_TP(type,1);
    sumTP_45 = @(type) w_term_TP(type,2) + w_term_TP(type,3) ;
    metric = @(c) (w_term_1(c,0)-w_term_TP(c,0))./w_term_TP(c,0);
    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{TP}})/{ \log \hat{W}_{TP}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-1 1]);
    exportgraphics(f,"../../Output/relative_12345_GFT_TP.pdf")
    mean(abs(metric(all)))

    metric = @(c) (w_term_2(c,0)-w_term_TP(c,0))./w_term_TP(c,0);
    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Constant} - \log \hat{W}_{TP}})/{ \log \hat{W}_{TP}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-1 1]);
    exportgraphics(f,"../../Output/crelative_12345_GFT_TP.pdf")

    metric = @(c) (w_term_1(c,0)-w_term_TP(c,0))./w_term_1(c,0);
    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{TP}})/{ \log \hat{W}_{Semi}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-1 1]);
    exportgraphics(f,"../../Output/crelative_12345_GFT_TP_1.pdf")
    
    metric = @(c) (w_term_1(c,0)-w_term_TP(c,0))./w_term_2(c,0);
    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{TP}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-1.2 1]);
    exportgraphics(f,"../../Output/crelative_12345_GFT_TP_2.pdf")
    

%% From the Literature (Log-Normal)
% Generate Log-Normal Data
% Fernandes, A., Klenow, P., Denisse Periola, M., Meleshchuk, S., and Rodriguez-Clare, A. (2017). The intensive margin in trade: Moving beyond pareto. Unpublished manuscript.
% https://siepr.stanford.edu/sites/default/files/publications/555wp_1.pdf
% Meyer et al
% https://pubs.aeaweb.org/doi/pdfplus/10.1257/aer.104.5.310

    n_sigma         = .797; % Table 1 Col 1 (From AER P&P) DO WITH 0.6 too
    % n_sigma         = .6; % Table 1 Col 1 (From AER P&P) DO WITH 0.6 too
    n_mu            = 0;
    % sigma_LN_Lit    = sigma;

    % original
    z_LN             = @(n) exp(n_mu + sqrt(2)*n_sigma*erfinv(2*(1-n)-1))+1;
    rho_bar_LN       = @(n)  (1/n).*integral(@(x) B.*(tau./(1+exp(n_mu + sqrt(2)*n_sigma.*erfinv(2*(1-x)-1)))).^(1-sigma) ,0,(n));
    rho_LN           = @(n) B.*(tau./z_LN(n)).^(1-sigma);
    epsilon_LN       = @(n) rho_LN(n);
    ratio_LN         = @(n) rho_LN(n)./epsilon_LN(n);
    C                = size(n_ij,1);
    top1             = arrayfun(@(lb,ub) integral(ratio_LN,lb,ub),zeros(C),n_ij);
    bot1             = arrayfun(@(lb,ub) integral(rho_LN,lb,ub)  ,zeros(C),n_ij)./arrayfun(epsilon_LN,n_ij);
    precalculated_integrals = top1./bot1; 
    FGFT = @(n_ii_hat) GFT(n_ii_hat,X_ij,epsilon_LN,rho_LN,rho_bar_LN,ratio_LN,n_ij,x_ij,1,precalculated_integrals,'array',1,0);
    [n_ii_hat_LN]    = lsqnonlin(FGFT,ones(size(n_ij,1),1)+.01,zeros(size(n_ij,1),1), 1./diag(n_ij),optimoptions('lsqnonlin','Display','iter'));
    [~,N_hat_LN]     = FGFT(n_ii_hat_LN);
    elasticity       = -1/(sigma-1);
    W_hat_splineLN_1 = diag(x_ij).^elasticity;
    W_hat_splineLN_2 = (n_ii_hat_LN.*N_hat_LN).^elasticity;
    W_hat_splineLN_3 = ( arrayfun(rho_bar_LN,n_ii_hat_LN.*diag(n_ij))./arrayfun(rho_bar_LN,diag(n_ij))    ).^elasticity;
    W_hat_splineLN   = W_hat_splineLN_1 .* W_hat_splineLN_2 .* W_hat_splineLN_3;
    
    LN_check = mean(abs(FGFT(n_ii_hat_LN)))
    % save('ESTIMATES_rewrite','W*','x_ij');

    results_LN = [log(W_hat_splineLN) log(W_hat_splineLN_1) log(W_hat_splineLN_2) log(W_hat_splineLN_3)];

    %% Plot LN
    w_term_LN = @(flag,i) results_LN(flag,(i+1)) ;
    sumLN_13 = @(type) w_term_LN(type,1);
    sumLN_45 = @(type) w_term_LN(type,2) + w_term_LN(type,3) ;
    metric = @(c) (w_term_1(c,0)-w_term_LN(c,0))./w_term_LN(c,0);
    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{LN}})/{ \log \hat{W}_{LN}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-1 1]);
    exportgraphics(f,"../../Output/relative_12345_GFT_LN.pdf")
    % exportgraphics(f,"../../Output/relative_12345_GFT_LN06.pdf")
    mean(abs(metric(all)))

    share = mean(n_ij-diag(diag(n_ij)),2);
    ln_n_ij = log(share);

    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Constant} - \log \hat{W}_{LN}})/{ \log \hat{W}_{LN}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-1 1]);
    exportgraphics(f,"../../Output/crelative_12345_GFT_LN.pdf")

    metric = @(c) (w_term_1(c,0)-w_term_LN(c,0))./w_term_1(c,0);
    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{LN}})/{ \log \hat{W}_{Semi}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-2 1]);
    exportgraphics(f,"../../Output/relative_12345_GFT_LN_1.pdf")

    metric = @(c) (w_term_1(c,0)-w_term_LN(c,0))./w_term_2(c,0);
    GFT_agg_plot
    ylabel({'$({\log \hat{W}_{Semi} - \log \hat{W}_{LN}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');        
    ylim([-1.2 1]);
    exportgraphics(f,"../../Output/relative_12345_GFT_LN_2.pdf")

    %% Text Points
        rich_deviation = mean((metric(rich)))
        poor_deviation = mean((metric(poor)))
        min((metric(poor)))
        [i1,i2] = min((metric(poor)));ll = labels(poor);ll(i2)
        
        max((metric(rich)))
        [i1,i2] = max((metric(rich)));ll = labels(rich);ll(i2)
          

    % %% Baseline
    % mean_n_ij   = mean(n_ij_init-diag(diag(n_ij_init)),2);
    % ln_n_ij     = log(mean_n_ij);

    % %% Alternative n_ij (weighted)
    % w_ij = X_ij_init-diag(diag(X_ij_init));
    % mean_n_ij   = sum(w_ij.*n_ij_init-diag(diag(w_ij.*n_ij_init)),2)./sum(w_ij,2);
    % ln_n_ij     = log(mean_n_ij);

    scatter(w_term_LN(all,0),w_term_2(all,0))
    hold on
    plot(w_term_LN(all,0),w_term_LN(all,0))
    hold off
    xlabel("Log Normal Welfare Gain",'FontSize',16,'Interpreter','Latex')
    ylabel("Constant Elasticitiy Welfare Gain",'FontSize',16,'Interpreter','Latex')

    xmetric  = @(c) w_term_1(c,0)
    ymetric1 = @(c) w_term_LN(c,0);
    a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    % text(xmetric(rich),ymetric1(rich),labels(rich),'VerticalAlignment','top','HorizontalAlignment','left')
    b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
    % text(xmetric(poor),ymetric1(poor),labels(poor),'VerticalAlignment','bottom','HorizontalAlignment','right')
    plot(xmetric(all),xmetric(all))
    hold off
    legend([a b ],{'Developed', 'Developing' },'FontSize',14,'Location','southeast')
    ylabel("Log Normal Welfare Gain $\hat{W}$",'FontSize',16,'Interpreter','Latex')
    xlabel("Semiparametric Welfare Gain $\hat{W}$",'FontSize',16,'Interpreter','Latex')
    exportgraphics(f,"../../Output/LN_SP.pdf")

    xmetric  = @(c) w_term_1(c,0)
    ymetric1 = @(c) w_term_TP(c,0);
    a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    % text(xmetric(rich),ymetric1(rich),labels(rich),'VerticalAlignment','top','HorizontalAlignment','left')
    b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
    % text(xmetric(poor),ymetric1(poor),labels(poor),'VerticalAlignment','bottom','HorizontalAlignment','right')
    plot(xmetric(all),xmetric(all))
    hold off
    legend([a b ],{'Developed', 'Developing' },'FontSize',14,'Location','southeast')
    ylabel("Truncated Pareto Welfare Gain $\hat{W}$",'FontSize',16,'Interpreter','Latex')
    xlabel("Semiparametric Welfare Gain $\hat{W}$",'FontSize',16,'Interpreter','Latex')
    exportgraphics(f,"../../Output/TP_SP.pdf")

    xmetric  = @(c) w_term_2(c,0)
    ymetric1 = @(c) w_term_LN(c,0);
    a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    % text(xmetric(rich),ymetric1(rich),labels(rich),'VerticalAlignment','top','HorizontalAlignment','left')
    b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
    % text(xmetric(poor),ymetric1(poor),labels(poor),'VerticalAlignment','bottom','HorizontalAlignment','right')
    plot(xmetric(all),xmetric(all))
    hold off
    legend([a b ],{'Developed', 'Developing' },'FontSize',14,'Location','southeast')
    ylabel("Log-Normal Welfare Gain $\hat{W}$",'FontSize',16,'Interpreter','Latex')
    xlabel("Constant Elasticity Welfare Gain $\hat{W}$",'FontSize',16,'Interpreter','Latex')
    f = gcf;
    exportgraphics(f,"../../Output/LN_CE.pdf")

%% Wrapper to Plot Aggregate Data
    xmetric  = @(c) w_term_1(c,0)
    ymetric1 = @(c) w_term_2(c,0);

    all = poor | rich;   
    a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    % text(xmetric(rich),ymetric1(rich),labels(rich),'VerticalAlignment','top','HorizontalAlignment','left')
    b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
    % text(xmetric(poor),ymetric1(poor),labels(poor),'VerticalAlignment','bottom','HorizontalAlignment','right')
    plot(xmetric(all),xmetric(all))
    hold off
    legend([a b ],{'Developed', 'Developing' },'FontSize',14,'Location','southeast')
    f = gcf;
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('Semiparametric $\hat{W}$','FontSize',16,'Interpreter','Latex');            
    ylabel('Constant Elasticity $\hat{W}$','FontSize',16,'Interpreter','Latex')
    f =gca;
    f.XTick   = [.001 .005 .01 .02 .05 .1 .2] ; f.XTickLabel = {'0.1%', '.5%' ,'1%', '2%','5%' '10%' '20%'};
    f.YTick   = [.001 .005 .01 .02 .05 .1 .2] ; f.YTickLabel = {'0.1%', '.5%' ,'1%', '2%','5%' '10%' '20%'};
    xlim([.007 .24]);        ylim([.007 .24]);
    exportgraphics(f,"../../Output/GFT_COMP.pdf")


    a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    % text(xmetric(rich),ymetric1(rich),labels(rich),'VerticalAlignment','top','HorizontalAlignment','left')
    b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
    % text(xmetric(poor),ymetric1(poor),labels(poor),'VerticalAlignment','bottom','HorizontalAlignment','right')
    plot(xmetric(all),xmetric(all))
    hold off
    legend([a b ],{'Developed', 'Developing' },'FontSize',14,'Location','southeast')
    xlabel('Semiparametric $\hat{W}$','FontSize',16,'Interpreter','Latex');            
    ylabel('Constant Elasticity $\hat{W}$','FontSize',16,'Interpreter','Latex')
    f =gca;
    exportgraphics(f,"../../Output/GFT_COMP_l.pdf")


     GFT_Script_v17_EKK

    