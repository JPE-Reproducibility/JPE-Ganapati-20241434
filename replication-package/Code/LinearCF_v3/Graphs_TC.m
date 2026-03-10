


%% Decomposition Exact - Elasticity
    dle = log(1+results_splin6.shock);
    fprintf('Spline 4-Group     \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',100*sum(results_splin6.Decomp_5.*weights)/dle)
    fprintf('  Developed Origins\t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',100*sum(diag(results_splin4.G_ij).*results_splin6.Decomp_5.*weights)/sum(diag(results_splin4.G_ij.*weights))/dle)
    fprintf('  Developing Origins\t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',100*sum(diag(1-results_splin4.G_ij).*results_splin6.Decomp_5.*weights)/sum(weights.*diag(1-results_splin4.G_ij))/dle)

    filename = "../../Output/Decomposition_Results_TC" + scenario+ ".csv";
    ROOT = "../..";
    fileID = fopen(filename, 'w');
    fprintf(fileID, 'Spline 4-Group     \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',    100 * sum(results_splin6.Decomp_5 .* weights) / dle);
    fprintf(fileID, '  Developed Origins\t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',    100 * sum(diag(results_splin4.G_ij) .* results_splin6.Decomp_5 .* weights) / sum(diag(results_splin4.G_ij .* weights)) / dle);
    fprintf(fileID, '  Developing Origins\t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',    100 * sum(diag(1 - results_splin4.G_ij) .* results_splin6.Decomp_5 .* weights) / sum(weights .* diag(1 - results_splin4.G_ij)) / dle);
    fclose(fileID);
    disp(['Output written to file: ' filename]);

%% Baseline
    mean_n_ij   = mean(n_ij_init-diag(diag(n_ij_init)),2);
    ln_n_ij     = log(mean_n_ij);

    % %% Alternative n_ij (weighted)
    % w_ij = X_ij_init-diag(diag(X_ij_init));
    % mean_n_ij   = sum(w_ij.*n_ij_init-diag(diag(w_ij.*n_ij_init)),2)./sum(w_ij,2);
    % ln_n_ij     = log(mean_n_ij);

%% Metrics
%     % WORLD_WELFARE = sum(results_linear.Decomp_5(:,1).*weights);
    % w_term = @(flag,i) results_splin6.Decomp_5(flag,(i+1)) - results_linear.Decomp_5(flag,(i+1)) ;
    IND = results_splin6.G_ij+eye(size(results_splin6.G_ij))*4;
    C3 = [.5 .25  0];         
    C1 = [.5 .5 1];           
    C2 = max(min(1,C1-.1),0);
    C4 = max(min(1,C3-.1),0);

    n_ii_hat        = diag(results_splin6.dn_linear);
    N_i_hat         = results_splin6.dN_linear;
    all             = poor | rich;

%% Relative Outcomes
    type0 = 'Developed -> Developed';
    type1 = 'Developed -> Developing';
    type2 = 'Developing -> Developed';
    type3 = 'Developing -> Developing';
    
    n = log(results_linear.n_ij_init);
    x = log(results_linear.x_ij_init);
    
    X_ij_hat_linear = results_linear.X_ij./results_linear.X_ij_init;
    X_ij_hat_spline = results_splin6.X_ij./results_splin6.X_ij_init;
    
    lxbar_ij_ratio = log(results_splin6.dxbar_linear ./results_linear.dxbar_linear)./log(results_linear.dxbar_linear);
    ln_ij_ratio    = log(results_splin6.dn_linear ./results_linear.dn_linear)./log(results_linear.dn_linear);
    lX_ij_ratio    = log(X_ij_hat_spline./X_ij_hat_linear)./log(X_ij_hat_linear);

%% lxbar_ij_ratio
    plot_var = lxbar_ij_ratio;
    CF_flow_plot
    ylabel({'${\log {\hat{X}^{Semi}_{ij}}}/{\log {\hat{X}^{Constant}_{ij}}} - 1$'},'FontSize',16,'Interpreter','Latex')
    ylim(([-1 2]));
    exportgraphics(f,[ '../../Output/'  'X_' results_linear.file '.pdf' ])

%% lxbar_ij_ratio
    plot_var = lxbar_ij_ratio;
    CF_flow_plot
    ylabel({'${\log {\bar{\hat{x}_{ij}}^{Semi}_{ij}}}/{\log {\bar{\hat{x}_{ij}}^{Constant}_{ij}}} - 1$'},'FontSize',16,'Interpreter','Latex')
    ylim(([-1 1]));
    exportgraphics(f,['../../Output/'  'xbar_' results_linear.file '.pdf' ])

%% ln_ij_ratio
    plot_var = lX_ij_ratio;
    CF_flow_plot
    ylabel({'${\log {\hat{n}^{Semi}_{ij}}}/{\log {\hat{n}^{Constant}_{ij}}} - 1$'},'FontSize',16,'Interpreter','Latex')
    ylim(([-1 2]));
    exportgraphics(f,['../../Output/'  'n' results_linear.file '.pdf' ])

%% Excel Diagnostics
if scenario == 2

    excel_output = "CF_results_revision_" + scenario + ".xlsx";

    Results_main = [  diag(x_ij_init) mean_n_ij results_splin6.Decomp_5 diag(results_splin4.G_ij) results_splin6.dN_linear diag(results_splin6.dn_linear) ln_x_ii ln_n_ij n_ii_hat N_i_hat];
    T = array2table(Results_main,'VariableNames', {'x_ii' ,'mean n_ij', 'Welfare', ...
        'Term1 A_tech' 'Term2 B_tot', ...
        'Term3 C_demand','Term4 D_Extensive','Term5 E_Selection','Error', ...
        'Developed Status','hat N_i','hat n_ii','ln_x_ii','ln_n_ij','n_ii_hat','N_i_hat'});
    T = [labels.textdata T];
    writetable(T,excel_output,'sheet','Spline OD');

    Results_main = [  diag(x_ij_init) mean_n_ij results_linear.Decomp_5 diag(results_splin4.G_ij) results_linear.dN_linear diag(results_linear.dn_linear) ln_x_ii ln_n_ij  ];
    T = array2table(Results_main,'VariableNames', {'x_ii' ,'mean n_ij', 'Welfare', ...
        'Term1 A_tech' 'Term2 B_tot', ...
        'Term3 C_demand','Term4 D_Extensive','Term5 E_Selection','Error', ...
        'Developed Status','hat N_i','hat n_ii','ln_x_ii','ln_n_ij'});
    T = [labels.textdata T];
    writetable(T,excel_output,'sheet','Constant Elasticitiy');

    %% Terms
        w_term_1 = @(flag,i) results_splin6.Decomp_5(flag,(i+1)) ;
        w_term_2 = @(flag,i) results_linear.Decomp_5(flag,(i+1)) ;
        sum1_13 = @(type) w_term_1(type,1) + w_term_1(type,2) + w_term_1(type,3);
        sum2_13 = @(type) w_term_2(type,1) + w_term_2(type,2) + w_term_2(type,3);
        sum1_45 = @(type) w_term_1(type,4) + w_term_1(type,5) ;
        sum2_45 = @(type) w_term_2(type,4) + w_term_2(type,5) ;
    

    %% Correlations
        A1 = corr(ln_n_ij(rich),(w_term_1(rich,0)-w_term_2(rich,0))./w_term_2(rich,0))
        A2 = corr(ln_n_ij(poor),(w_term_1(poor,0)-w_term_2(poor,0))./w_term_2(poor,0))
    

    %% Text Figures
        metric = @(c) (w_term_1(c,0)-w_term_2(c,0))./w_term_2(c,0);
        mean(metric(rich))
        mean(metric(poor))

        rich_deviation = mean(abs(metric(rich)))
        poor_deviation = mean(abs(metric(poor)))

        [i1,i2] = sort((metric(rich)))
        % ll(i2)

        [i1,i2] = sort((metric(poor)))
        % ll(i2)

        
    %% Output Plots For 1-2-3-4-5
        metric = @(c) (w_term_1(c,0)-w_term_2(c,0))./w_term_2(c,0);
        CF_agg_plot
        ylabel({'$(\log \hat{W}_{Semi} - \log \hat{W}_{Constant})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
        legend([a b ],{'Developed Countries', 'Developing Countries'},'FontSize',14,'Location','best')
        exportgraphics(f, ROOT + "/Output/trelative_12345_" + results_linear.file + ".pdf")
    
    %% Output Plots For 1-2-3
        metric = @(c) (sum1_13(c)-sum2_13(c))./w_term_2(c,0);
        CF_agg_plot
        legend([a b ],{'Developed Countries', 'Developing Countries'},'FontSize',14,'Location','best')
        ylabel({'$({\log \hat{W}^{Neoclassical}_{Semi} - \log \hat{W}^{Neoclassical}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/trelative_123_" + results_linear.file + ".pdf")
    
    %% Output Plots For4-5
        metric = @(c) (sum1_45(c)-sum2_45(c))./w_term_2(c,0);
        CF_agg_plot
        legend([a b ],{'Developed Countries', 'Developing Countries'},'FontSize',14,'Location','best')
        ylabel({'$({\log \hat{W}^{Firm}_{Semi} - \log \hat{W}^{Firm}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/trelative_45_" + results_linear.file + ".pdf")
    
    %% Diagnostic Appendix Plots
    
    %% Bar Chart of Welfare Changes
        barh(results_splin6.Decomp_5(all,2:6)./results_splin6.Decomp_5(all,1), 'stacked')
        legend({'Technology', 'Terms of Trade','Demand Substitution','Firm Entry','Firm Selection'},'FontSize',12,'Location','northwest')
        yticklabels(labels.textdata(all))
        yticks(1:85)
        file =  ROOT + "/Output/Decomposition" + results_linear.file + ".pdf";
        f =gca;
        set(gcf,'position',[10,10,800,1000])
        exportgraphics(f,file)
        clf
        close all



    %% Scatter N_i
        plot_var = log(results_splin6.dN_linear)./log(1+results_splin6.shock);
        CF_diagnostic
        % ylabel({'$\hat{N}_i$'},'FontSize',16,'Interpreter','Latex');           
        % ylabel({'$Elasticity: \partial \hat{N}_i / \partial {shock}$'},'FontSize',16,'Interpreter','Latex');           
        ylabel({'$Elasticity:  {N}_i $'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/Ni_Hat" + results_linear.file + ".pdf")
    


        
    %% Scatter n_ii
        plot_var = log(n_ii_hat)./log(1+results_splin6.shock);
        CF_diagnostic
        ylabel({'$Elasticity: {n}_{ii}$'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/nii_Hat" + results_linear.file + ".pdf")
    
    
    %% Scatter n_ij
        w_ij = X_ij_init-diag(diag(X_ij_init));
        ln_ij_hat_w = w_ij.*log(results_splin6.dn_linear);
        ln_ij_hat        = sum(ln_ij_hat_w,2)./sum(w_ij,2);
        plot_var = ln_ij_hat./log(1+results_splin6.shock);
        CF_diagnostic
        ylabel({'$Mean_{j \ne i}^{weight=x_{ij}^{0}} \log  \hat{n}_{ij} / \log {shock}$'},'FontSize',16,'Interpreter','Latex');           
        ylabel({'Mean Elasticity$ _{j \ne i}^{weight=x_{ij}^{0}}   {n}_{ij} $'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/nij_Hat" + results_linear.file + ".pdf")
    

    
    %% Scatter n_ji
        w_ij = X_ij_init-diag(diag(X_ij_init));
        ln_ij_hat_w = w_ij.*log(results_splin6.dn_linear);
        ln_ji_hat        = sum(ln_ij_hat_w,1)./sum(w_ij,1);
        plot_var = ln_ji_hat./log(1+results_splin6.shock);
        CF_diagnostic
        ylabel({'$Mean_{i \ne j}^{weight=x_{ij}^{0}} \log  \hat{n}_{ij} / \log {shock}$'},'FontSize',16,'Interpreter','Latex');           
        ylabel({'Mean Elasticity$ _{i \ne j}^{weight=x_{ij}^{0}}   {n}_{ij} $'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/nji_Hat" + results_linear.file + ".pdf")



    %% Scatter n_ji
        w_ij = X_ij_init;
        ln_ij_hat_w = w_ij.*log(results_splin6.dn_linear);
        ln_ji_hat        = sum(ln_ij_hat_w,1)./sum(w_ij,1);
        plot_var = ln_ji_hat./log(1+results_splin6.shock);
        CF_diagnostic
        ylabel({'$Mean_{i \ne j}^{weight=x_{ij}^{0}} \log  \hat{n}_{ij} / \log {shock}$'},'FontSize',16,'Interpreter','Latex');           
        ylabel({'Mean Elasticity$ _{i}^{weight=x_{ij}^{0}}   {n}_{ij} $'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/nji2_Hat" + results_linear.file + ".pdf")

    %% Scatter n_ji
        w_ij = X_ij_init;
%         ln_ij_hat_w = w_ij.*log(results_splin6.dn_linear);
%         ln_ji_hat        = sum(ln_ij_hat_w,1)./sum(w_ij,1);
        plot_var = results_splin6.Decomp_5(:,6)./log(1+results_splin6.shock);
        CF_diagnostic
        ylabel({'$Mean_{i \ne j}^{weight=x_{ij}^{0}} (\rho_{ij}^0 + 1 )\log  \hat{n}_{ij} / \log {shock}$'},'FontSize',16,'Interpreter','Latex');           
        ylabel({'Mean Elasticity$ _{i \ne j}^{weight=x_{ij}^{0} (\rho_{ij}^0 + 1 )}   {n}_{ij} $'},'FontSize',16,'Interpreter','Latex');           
        exportgraphics(f,ROOT + "/Output/rhoji_Hat" + results_linear.file + ".pdf")



end