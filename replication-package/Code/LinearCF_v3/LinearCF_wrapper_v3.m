%% Run Piece-wise Linear
    clear;
    close all
    addpath( '../GMM_estimation_v3')

%% Run Simulation for Counterfactuals
% CF are currently hard-coded
    Segments = 1100;
    sigma    = 3.2;

    for scenario = [ 5 2]
        
        [Real_Wage_linear,dN_linear,~,~,results_linear] = LinearCF( Segments,sigma,scenario,0);
        [Real_Wage_splin4,dN_splin4,~,~,results_splin4] = LinearCF( 2,sigma,scenario,4);
        [Real_Wage_splin6,dN_splin6,labels,~,results_splin6] = LinearCF( Segments,sigma,scenario,6);

        %% Setup Output Files
        n_ij_init   = results_splin6.n_ij_init;
        x_ij_init   = results_linear.x_ij_init;
        X_ij_init   = results_linear.X_ij_init;
        weights     = sum(results_linear.X_ij_init,1)';
        weights     = weights/sum(weights);
        x_ii        = diag(x_ij_init);
        mean_n_ij   = mean(n_ij_init-diag(diag(n_ij_init)),2);
        ln_x_ii     = log(1./(diag(x_ij_init)));
        ln_n_ij     = log(mean_n_ij);

        sample_selection =    csvread("../../Data/Int/WIOD_sampleB/Sample_2012B.csv");
        
        OD = diag(results_splin6.G_ij);
        rich = (OD==0) & (sample_selection == 3);
        poor = (OD==3) & (sample_selection == 3);

        group1 = rich;
        group2 = poor;

        if scenario == 4 || scenario == 5
            Graphs_TC
            Graphs_GSP
        end

        if scenario == 1 || scenario == 2 || scenario == 3 || scenario == 999
            Graphs_TC
        end

    end
