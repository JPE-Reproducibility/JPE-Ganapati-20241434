%% Setup
    GSP_donor = sum(results_splin6.GSP_ij,1)'>10 & (sample_selection == 3) ;
    GSP_recip = sum(results_splin6.GSP_ij,2)>10 & GSP_donor == 0 & (sample_selection == 3);
    GSP_other = GSP_donor == 0 & GSP_recip == 0 & (sample_selection == 3);

    group1 = GSP_donor;
    group2 = GSP_recip;
            
    w_term_1 = @(flag,i) results_splin6.Decomp_5(flag,(i+1)) ;
    w_term_2 = @(flag,i) results_linear.Decomp_5(flag,(i+1)) ;
    sum1_13 = @(type) w_term_1(type,1) + w_term_1(type,2) + w_term_1(type,3);
    sum2_13 = @(type) w_term_2(type,1) + w_term_2(type,2) + w_term_2(type,3);
    sum1_45 = @(type) w_term_1(type,4) + w_term_1(type,5) ;
    sum2_45 = @(type) w_term_2(type,4) + w_term_2(type,5) ;

%% Decomposition Elasticity - Weighted
    dle = log(1+.01);
    fprintf('\n\nDecomposition Exact - UW\n')
    fprintf('Spline 4-Group     \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',100*sum(results_splin6.Decomp_5.*weights)/sum(weights)/dle)
    fprintf('  GSP Donors     \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',100*sum(GSP_donor.*results_splin6.Decomp_5.*weights)/sum(GSP_donor.*weights)/dle)
    fprintf('  GSP Recipients \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',100*sum(GSP_recip.*results_splin6.Decomp_5.*weights)/sum(GSP_recip.*weights)/dle)
    fprintf('  Non-GSP        \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',100*sum(GSP_other.*results_splin6.Decomp_5.*weights)/sum(GSP_other.*weights/dle))

    filename = "../../Output/Decomposition_Results_GSP" + scenario+ ".csv";
    fileID = fopen(filename, 'w');
    fprintf(fileID, '\n\nDecomposition Exact - UW\n');
    fprintf(fileID, 'Spline 4-Group     \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',     100 * sum(results_splin6.Decomp_5 .* weights) / sum(weights) / dle);
    fprintf(fileID, '  GSP Donors     \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',     100 * sum(GSP_donor .* results_splin6.Decomp_5 .* weights) / sum(GSP_donor .* weights) / dle);
    fprintf(fileID, '  GSP Recipients \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',     100 * sum(GSP_recip .* results_splin6.Decomp_5 .* weights) / sum(GSP_recip .* weights) / dle);
    fprintf(fileID, '  Non-GSP        \t %10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\t%10.8f\n',     100 * sum(GSP_other .* results_splin6.Decomp_5 .* weights) / sum(GSP_other .* weights / dle));
    fclose(fileID);
    disp(['Output written to file: ' filename]);

%% Output Results In Excel
    excel_output = "CF_results_revision_" + scenario + ".xlsx";
    Results_main = [  diag(x_ij_init) mean_n_ij results_linear.Decomp_5 diag(results_splin4.G_ij) GSP_donor GSP_recip GSP_other weights];
    T = array2table(Results_main,'VariableNames', {'x_ii' ,'mean n_ij', 'Welfare', ...
        'Term1 A_tech' 'Term2 B_tot', ...
        'Term3 C_demand','Term4 D_Extensive','Term5 E_Selection','Error', ...
        'Developed Status','GSP Donor','GSP Recipient','GSP Other','weights'});
    T = [labels.textdata T];
    writetable(T,excel_output,'sheet','Generalized Pareto');
    Results_main = [  diag(x_ij_init) mean_n_ij results_splin6.Decomp_5 diag(results_splin4.G_ij) GSP_donor GSP_recip GSP_other];
    T = array2table(Results_main,'VariableNames', {'x_ii' ,'mean n_ij', 'Welfare', ...
        'Term1 A_tech' 'Term2 B_tot', ...
        'Term3 C_demand','Term4 D_Extensive','Term5 E_Selection','Error', ...
        'Developed Status','GSP Donor','GSP Recipient','GSP Other'});
    T = [labels.textdata T];
    writetable(T,excel_output,'sheet','Spline OD');

%% Output Plots For 1-2-3-4-5
    all = GSP_donor | GSP_recip;
    metric = @(c) (w_term_1(c,0)-w_term_2(c,0))./w_term_2(c,0);
    CF_agg_plot
    ylabel({'$({\log \hat{W}_{semi} - \log \hat{W}_{Constant}}/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
    legend([a b],{'Donors', 'Beneficiaries'},'FontSize',14,'Location','best')
    exportgraphics(f,ROOT+"/Output/relative_12345_" + results_linear.file + ".pdf")

    ylim([-3 1]);
    sort(metric(all))

%% Output Plots For 1-2-3
    metric =  @(c) (sum1_13(c)-sum2_13(c))./w_term_2(c,0);
    CF_agg_plot
    ylabel({'$({\log \hat{W}^{Neoclassical}_{Semi} - \log \hat{W}^{Neoclassical}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
    legend([a b],{'Donors', 'Beneficiaries'},'FontSize',14,'Location','best')
    exportgraphics(f,ROOT + "/Output/relative_123_" + results_linear.file + ".pdf")

%% Output Plots For 45
    metric =  @(c) (sum1_45(c)-sum2_45(c))./w_term_2(c,0);
    CF_agg_plot
    ylabel({'$({\log \hat{W}^{Firm}_{Semi} - \log \hat{W}^{Firm}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
    legend([a b],{'Donors', 'Beneficiaries'},'FontSize',14,'Location','best')
    exportgraphics(f,ROOT + "/Output/relative_45_" + results_linear.file + ".pdf")

    ylim([-5 5]);
    sort(metric(all))

%% Output Plots For 4
    metric =  @(c) (w_term_1(c,4)-w_term_2(c,4))./w_term_2(c,0);
    CF_agg_plot
    ylabel({'$({\log \hat{W}^{Firm}_{Entry} - \log \hat{W}^{Firm}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
    legend([a b],{'Donors', 'Beneficiaries'},'FontSize',14,'Location','best')
    ylim([-0.5 2.5]);
    exportgraphics(f,ROOT + "/Output/relative_4_" + results_linear.file + ".pdf")

%% Output Plots For 5
    all = GSP_donor | GSP_recip;
    metric =  @(c) (w_term_1(c,5)-w_term_2(c,5))./w_term_2(c,0);
    CF_agg_plot
    ylabel({'$({\log \hat{W}^{Firm}_{Selection} - \log \hat{W}^{Firm}_{Constant}})/{ \log \hat{W}_{Constant}}$'},'FontSize',16,'Interpreter','Latex');           
    legend([a b],{'Donors', 'Beneficiaries'},'FontSize',14,'Location','best')
    ylim([-2 1]);
    exportgraphics(f,ROOT + "/Output/relative_5_" + results_linear.file + ".pdf")


