%% Set Directories
    clear
    close all
    data_dir = '../../Data/Int/WIOD_sampleB/';
    disp(data_dir)
    mkdir('output')
    
    %% Set up parameter space
    x               = -10:.01:-.01;
    C1 = [1 0  0];            C2 = [.5 0 0];
    C3 = [0 .5 1];            C4 = [0 .5 .5];
    z(1).color = [.5 .5 .5];      z(2).color = [.5 .75 .5];
    z(3).color = [.75 .5 .15];    z(4).color = [.75 .15 .25];
    types = {'rho','epsilon','theta','intensive','extensive'};

    kappa_tau = 1;
    kappa_f   = 0;
    sigma     = 3.2;

    kappa_epsilon = 1/((sigma-1)*kappa_tau+kappa_f);
    gamma_guess = [.2 0 0 0 0 .2 0 0 0 0];

    basefile = '2012_Tetiuw.csv'; name = '2012_Teti';  grav_spec = 'all';

%% ADD EXPLANATORY VARIABLES
    fix = 2;
    knots = 3;
    M_base = csvread(strcat(data_dir,basefile));
    d = make_data_combo(M_base,knots,'cubic','','',0,grav_spec);    d.sigma = sigma;       d.x = x;
    [~,o] =  GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,'base',fix,gamma_guess,grav_spec );
    start = size(o.ests2,2)-2*size(o.est_epsilon,2)-1;
    gamma_gravity = o.ests2((end-start):end)
    save 'gamma_gravity' gamma_gravity
    p_overlay_elasticity( d,o,'theta'       ,'../../Output/full_SE' );

    o.title = '+ Gravity (FTA + common curr + language + colony)';     

%% Run Programs
    Simulation_Data_v10
    GenerateBaseline

%% Baseline (taking the first stage results as given)
    % Take Pass-through as a given
    [~,o] =  GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,'base',0,gamma_gravity,grav_spec );
    Figure_Sequencer(d,o,'F2a_base','')
    o.title = '+ Gravity (FTA + common curr + language + colony)';     

    figure;
    p_overlay_elasticity( d,o,'epsilon_CC'   ,'../../Output/F2a_base_' );
    p_overlay_elasticity( d,o,'rho_CC'       ,'../../Output/F2a_base_' );
    p_overlay_elasticity( d,o,'extensive_CC'       ,'../../Output/F2a_base_' );
    h01 = p_overlay_elasticity( d,o, 'theta', '../../Output/F2a_base_'   );



%% IV Testing
    M_IV = csvread(strcat(data_dir,'2012_TuwIV2.csv'));
    dIV = make_data_combo(M_IV,knots,'cubic','IV','',0,'dist');  dIV.sigma = sigma;  dIV.x = x; 
    [~,oIV] =  GMM_wrapper_gravity( dIV, kappa_tau, kappa_epsilon,'base',2,[1 1],'dist' );
    p_overlay_elasticity( dIV,oIV,'theta'       ,'../../Output/F5_IV_' );
    start = size(oIV.ests2,2)-2*size(oIV.est_epsilon,2)-1;
    gamma_gravityIV = oIV.ests2((end-start):end);

    o.title = 'Baseline';     
    oIV.title = 'Instrumental Variables';      oIV.standard_errors1 = []; oIV.standard_errors2 = [];
    types2 = {'rho','epsilon','theta','intensive','extensive'};
    for t = 1:1:size(types2,2)
        graph_overlay([ o oIV],[ d dIV],types2{t},'../../Output/F5_IV_',  "", 0)
    end

%%  1 knot
    fix = 0;
    knots = 1;
    d2q_Wo1 = make_data_combo(M_base,knots,'cubic','','',0,'all'); 
    d2q_Wo1.sigma = sigma;  d2q_Wo1.x = x;
    [~,o2q_Wo1] =  GMM_wrapper_gravity( d2q_Wo1, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    Figure_Sequencer(d2q_Wo1,o2q_Wo1,'F1_pooled','')

%%  1 knot Split
    fix = 0;
    knots = 1;
    d2q_Wo1 = make_data_combo(M_base,knots,'cubic','wealth_oI','',0,'all'); 
    d2q_Wo1.sigma = sigma;  d2q_Wo1.x = x;
    [~,o2q_Wo1] =  GMM_wrapper_gravity( d2q_Wo1, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    Figure_Sequencer(d2q_Wo1,o2q_Wo1,'F3a_origin_one','_split')


%% Wealth Origin
    fix = 0;
    % knots = 34;
    knots = 33;
    d2q_Wo = make_data_combo(M_base,knots,'cubic','wealth_oI','',0,'all'); 
    d2q_Wo.sigma = sigma;  d2q_Wo.x = x;
    [~,o2q_Wo] =  GMM_wrapper_gravity( d2q_Wo, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    p_overlay_elasticity( d2q_Wo,o2q_Wo,'theta_split'       ,'../../Output/test' );
    Figure_Sequencer(d2q_Wo,o2q_Wo,'F3a_origin','_split')


%% Wealth Origin with IV
    dIVs = make_data_combo(M_IV,knots,'cubic','wealth_oIV','',0,'dist');  dIVs.sigma = sigma;  dIVs.x = x; 
    [~,oIVs] =  GMM_wrapper_gravity( dIVs, kappa_tau, kappa_epsilon,'base',0,gamma_gravityIV,'dist' );
    p_overlay_elasticity( dIVs,oIVs,'theta_split'       ,'../../Output/F5_IVo_' );
    
        


%% Wealth Destination
% Simultaneous
    fix = 0;
    knots = 33;
    d2q_Wd = make_data_combo(M_base,knots,'cubic','wealth_dI','',0,'all'); 
    d2q_Wd.sigma = sigma;  d2q_Wd.x = x;
    [~,o2q_Wd] =  GMM_wrapper_gravity( d2q_Wd, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    Figure_Sequencer(d2q_Wd,o2q_Wd,'F3d2_destination','_split')

%% Deep Integration
% Simultaneous
% Something is off here
    fix = 0;
    knots = 33;
    d2q_comcur = make_data_combo(M_base,knots,'cubic','deep','',0,'all'); 
    d2q_comcur.sigma = sigma;  d2q_comcur.x = x;
    [~,o2q_comcur] =  GMM_wrapper_gravity( d2q_comcur, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    Figure_Sequencer(d2q_comcur,o2q_comcur,'F3d4_Deep','_split')

%% Language or Colony
% Simultaneous
    fix = 0;
    knots = 33;
    d2q_LC = make_data_combo(M_base,knots,'cubic','langcol','',0,'all'); 
    d2q_LC.sigma = sigma;  d2q_LC.x = x;
    [~,o2q_LC] =  GMM_wrapper_gravity( d2q_LC, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    Figure_Sequencer(d2q_LC,o2q_LC,'F3d5_LC','_split')

%% Export estimation sample under the assumption of n_ii = 1
    M_alt = csvread(strcat(data_dir,'2012_survival_all_.csv'));
    knots = 3;
    fix = 0;
    base_alt = make_data_combo(M_alt,knots,'cubic','','',0,'all');     base_alt.sigma = sigma;       base_alt.x = x;
    [~,obase_alt] =  GMM_wrapper_gravity( base_alt, kappa_tau, kappa_epsilon,'base_alt',fix,gamma_gravity,'all');
    o.title = 'Baseline';     
    obase_alt.title = 'Assume n_{ii} = 1';
    types = {'rho','theta','extensive'};
    for t = 1:1:size(types,2)
        graph_overlay([ o obase_alt],[ d base_alt],types{t},'../../Output/F5_nii1',  "", 0)
    end


%% Export estimation sample without n_iie
    M_alt = csvread(strcat(data_dir,'2012_nonii_.csv'));
    knots = 3;
    fix = 0;
    base_alt = make_data_combo(M_alt,knots,'cubic','','',0,'all');     base_alt.sigma = sigma;       base_alt.x = x;
    [~,obase_alt] =  GMM_wrapper_gravity( base_alt, kappa_tau, kappa_epsilon,'base_alt',fix,gamma_gravity,'all');
    obase_alt.title = 'Drop n_{ii}';     
    for t = 1:1:size(types,2)
        graph_overlay([ o obase_alt],[ d base_alt],types{t},'../../Output/F5_nonii',  "", 0)
    end

%% Export estimation sample without imputed survival rate
    M_alt = csvread(strcat(data_dir,'2012_survival_.csv'));
    knots = 3;
    fix = 0;
    base_alt = make_data_combo(M_alt,knots,'cubic','','',0,'all');     base_alt.sigma = sigma;       base_alt.x = x;
    [~,obase_alt] =  GMM_wrapper_gravity( base_alt, kappa_tau, kappa_epsilon,'base_alt',fix,gamma_gravity,'all');
    obase_alt.title = 'No Imputed Survival Rate';     
    for t = 1:1:size(types,2)
        graph_overlay([ o obase_alt],[ d base_alt],types{t},'../../Output/F5_noimpnii',  "", 0)
    end

%% Export estimation sample with 3 year rate
    M_alt = csvread(strcat(data_dir,'2012_survival_alt3_.csv'));
    knots = 3;      fix = 0;
    base_alt = make_data_combo(M_alt,knots,'cubic','','',0,'all');     base_alt.sigma = sigma;       base_alt.x = x;
    [~,obase_alt] =  GMM_wrapper_gravity( base_alt, kappa_tau, kappa_epsilon,'base_alt',fix,gamma_gravity,'all');
    obase_alt.title = '3-year survival rate';     
    for t = 1:1:size(types,2)
        graph_overlay([ o obase_alt],[ d base_alt],types{t},'../../Output/F5_nii3',  "", 0)
    end

%% Do Linear
    fix = 0;
    p_lin = make_data_combo(M_base,1,'cubic','','',0,'all');    p_lin.sigma = sigma;       p_lin.x = x;
    [~,o_lin] =  GMM_wrapper_gravity( p_lin, kappa_tau, kappa_epsilon,'base',fix,gamma_gravity,'all' );
    start = size(o_lin.ests2,2)-2*size(o_lin.est_epsilon,2)-1;
    gamma_linear = o_lin.ests2((end-start):end);
    h00 = p_overlay_elasticity( p_lin,o_lin, 'theta', '../../Output/test'   );

%% Use the Literature
    o.title = 'Baseline';     
    o_lin.title = 'Constant Elasticity';     
    for t = 1:1:size(types,2)
        graph_overlay([ o o_lin],[ d p_lin],types{t},'../../Output/F6_lit_',  "", 1)
    end


    

%% 4-way independent samples (this means different FE for each
% Wealth Origin-Destination
% Simultaneous
    fix = 0;
    knots = 3;
    d4q_Wod = make_data_combo(M_base,knots,'cubic','wealth_odI','',0,'all'); 
    d4q_Wod.sigma = sigma;  d4q_Wod.x = x;
    [~,o4q_Wod] =  GMM_wrapper_gravity( d4q_Wod, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    Figure_Sequencer(d4q_Wod,o4q_Wod,'F3d1_Cross','_quad')
    h03 = p_overlay_elasticity( d4q_Wod,o4q_Wod, 'theta', '../../Output/test'   );

    p_overlay_elasticity( d4q_Wod,o4q_Wod,  'theta_quad_nose'  ,'../../Output/F3d1_Cross_' );
    p_overlay_elasticity( d4q_Wod,o4q_Wod,  'theta_quadA'  ,'../../Output/F3d1_Cross_' );
    p_overlay_elasticity( d4q_Wod,o4q_Wod,  'theta_quadB'  ,'../../Output/F3d1_Cross_' );
    p_overlay_elasticity( d4q_Wod,o4q_Wod,  'extensive_quad'  ,'../../Output/F3d1_Cross_' );


%% 4-way independent samples (this means different FE for each
% pareto
    fix = 0;
    knots = 1;
    d4q_Wod_l = make_data_combo(M_base,knots,'cubic','wealth_odI','',0,'all'); 
    d4q_Wod_l.sigma = sigma;  d4q_Wod_l.x = x;
    [~,o4q_Wod_l] =  GMM_wrapper_gravity( d4q_Wod_l, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity,'all');
    h02 = p_overlay_elasticity( d4q_Wod_l,o4q_Wod_l, 'theta_quad', '../../Output/test'   );



%% 4-way independent samples (this means different FE for each
% Wealth Origin =- 4 categories
% Simultaneous
% gamma_alt_UW =    0.3161    0.2699
    fix = 0;
    knots = 3;
    d4oq_Wod = make_data_combo(M_base,knots,'cubic','wealth_ooI','',0,'all'); 
    d4oq_Wod.sigma = sigma;  d4oq_Wod.x = x;
    [~,oo4q_Wod] =  GMM_wrapper_gravity( d4oq_Wod, kappa_tau, kappa_epsilon,'test',fix,gamma_gravity,'all');
    Figure_Sequencer(d4oq_Wod,oo4q_Wod,'F3d5_4origin','_quad')


%% Graph distribution of results (Thetas)
    load('literature.mat')
    edges = 2:.1:9;
    aa = histogram(h00,edges,'Normalization','probability');
    hold on
    bb = histogram(h01,edges,'Normalization','probability');
    dd = histogram(h03,edges,'Normalization','probability');
    hold off
    legend([ aa bb  dd ],{'Constant Elasticity - 1 Group','Semiparametric - 1 Group', ...
        'Semiparametric - 4 Groups'},'FontSize',14,'NumColumns',1,'Location','northeast')
    xlabel("Distribution of \theta_{ij}",'FontSize',18);
    saveas(gcf,'../../Output/theta_dist_all2','epsc')




%% 4 knots
    M4 = csvread(strcat(data_dir,basefile));
    knots = 4; fix = 0;
    d_k4 = make_data_combo(M4,knots,'cubic','','',0,'all');     d_k4.sigma = sigma;       d_k4.x = x;
    [~,o_k4] =  GMM_wrapper_gravity( d_k4, kappa_tau, kappa_epsilon,'base_alt',fix,gamma_gravity,'all');
    o_k4.title = '4 knots';     
    for t = 1:1:size(types,2)
        graph_overlay([ o o_k4],[ d d_k4],types{t},'../../Output/F4_k4_',  "", 0)
    end



%% 5 knots
    M5 = csvread(strcat(data_dir,basefile));
    knots = 5; fix = 2;
    d_k5 = make_data_combo(M5,knots,'cubic','','',0,'all');     d_k5.sigma = sigma;       d_k5.x = x;
    [~,o_k5] =  GMM_wrapper_gravity( d_k5, kappa_tau, kappa_epsilon,'base_alt',fix,gamma_gravity,'all');
    o_k5.title = '5 knots';     
    for t = 1:1:size(types,2)
        graph_overlay([ o o_k5],[ d d_k5],types{t},'../../Output/F4_k5_',  "", 0)
    end



%% kappa_f = 1, gamma flexible
    knots = 3;
    fix = 2;
    kappa2_tau = 1;
    kappa2_f   = 1;
    kappa2_epsilon_f = 1/((sigma-1)*kappa2_tau+kappa2_f);
    basekf1 = make_data_combo(M_base,knots,'cubic','','',0,'dist');     basekf1.sigma = sigma;       basekf1.x = x;
    [~,obasekf1] =  GMM_wrapper_gravity( basekf1, kappa2_tau, kappa2_epsilon_f,'test',fix,gamma_gravity,'dist');
    obasekf1.title = '\kappa_e = 1';     
    for t = 1:1:size(types,2)
        graph_overlay([ o obasekf1],[ d basekf1],types{t},'../../Output/F5_basekf1',  "", 0)
    end
    start = size(obasekf1.ests2,2)-2*size(obasekf1.est_epsilon,2)-1;
    gamma_gravityf = obasekf1.ests2((end-start):end);


%% kappa_f = 1, gamma fixed - origin
    fix = 0;
    knots = 33;
    d2q_Wof = make_data_combo(M_base,knots,'cubic','wealth_oI','',0,'dist'); 
    d2q_Wof.sigma = sigma;  d2q_Wof.x = x;
    [~,o2q_Wof] =  GMM_wrapper_gravity( d2q_Wof, kappa2_tau, kappa2_epsilon_f,'test',fix ,gamma_gravityf,'dist');
    Figure_Sequencer(d2q_Wof,o2q_Wof,'F5_basekf1_origin','_split')


%% kappa_f = .5, gamma flexible
    knots = 3;
    fix = 2;
    kappa2_tau = 1;
    kappa2_f   = .5;
    kappa2_epsilon_f = 1/((sigma-1)*kappa2_tau+kappa2_f);
    basekf5 = make_data_combo(M_base,knots,'basekf5','','',0,'dist');     basekf5.sigma = sigma;       basekf5.x = x;
    [~,obasekf5] =  GMM_wrapper_gravity( basekf1, kappa2_tau, kappa2_epsilon_f,'test',fix,gamma_gravity,'dist');
    obasekf5.title = '\kappa_e = 1/2';     
    for t = 1:1:size(types,2)
        graph_overlay([ o obasekf5],[ d basekf5],types{t},'../../Output/F5_basekf5',  "", 0)
    end
    start = size(obasekf5.ests2,2)-2*size(obasekf5.est_epsilon,2)-1;
    gamma_gravityf = obasekf5.ests2((end-start):end);


    
%% kappa_f = .5, gamma fixed - origin
    fix = 0;
    knots = 33;
    d2q_Wof = make_data_combo(M_base,knots,'cubic','wealth_oI','',0,'dist'); 
    d2q_Wof.sigma = sigma;  d2q_Wof.x = x;
    [~,o2q_Wof] =  GMM_wrapper_gravity( d2q_Wof, kappa2_tau, kappa2_epsilon_f,'test',fix ,gamma_gravityf,'dist');
    Figure_Sequencer(d2q_Wof,o2q_Wof,'F5_basekf5_origin','_split')



%% Sigma 2.4
    fix = 2;
    knots = 3;
    ps_s24 = make_data_combo(M_base,knots,'cubic','','',0,'all');    ps_s24.sigma = 2.4;       ps_s24.x = x;
    [~,os_s24] =  GMM_wrapper_gravity( ps_s24, kappa_tau, kappa_epsilon,'base',fix,gamma_guess,'all' );
    start = size(os_s24.ests2,2)-2*size(os_s24.est_epsilon,2)-1;
    gamma_gravity_s24 = os_s24.ests2((end-start):end);

    % Baseline (taking the first stage results as given)
    % Take Pass-through as a given
    [~,os_s24] =  GMM_wrapper_gravity( ps_s24, kappa_tau, kappa_epsilon,'base',0,gamma_gravity_s24,'all' );
    os_s24.title = '\sigma = 2.4 (25th percentile)';     
    for t = 1:1:size(types,2)
        graph_overlay([ o os_s24],[ d ps_s24],types{t},'../../Output/F5_sigma24',  "", 0)
    end

%% Sigma 3.4
    fix = 2;
    knots = 3;
    ps_s34 = make_data_combo(M_base,knots,'cubic','','',0,'all');    ps_s34.sigma = 3.4;       ps_s34.x = x;
    [~,os_s34] =  GMM_wrapper_gravity( ps_s34, kappa_tau, kappa_epsilon,'base',fix,gamma_guess,'all' );
    start = size(os_s34.ests2,2)-2*size(os_s34.est_epsilon,2)-1;
    gamma_gravity_s34 = os_s34.ests2((end-start):end);


    % Baseline (taking the first stage results as given)
    % Take Pass-through as a given
    [~,os_s34] =  GMM_wrapper_gravity( ps_s34, kappa_tau, kappa_epsilon,'base',0,gamma_gravity_s34,'all' );
    os_s34.title = '\sigma = 3.4 (75th percentile)';     
    for t = 1:1:size(types,2)
        graph_overlay([ o os_s34],[ d ps_s34],types{t},'../../Output/F5_sigma34',  "", 0)
    end


%% 2010
    fix = 2;
    knots = 3;
    M2010 = csvread(strcat(data_dir,'2010_Tuw.csv'));
    type = 'all';
    ps_10 = make_data_combo(M2010,knots,'cubic','','',0,type);    ps_10.sigma = sigma;       ps_10.x = x;
    [~,os_10] =  GMM_wrapper_gravity( ps_10, kappa_tau, kappa_epsilon,'base',fix,gamma_guess,type );
    start = size(os_10.ests2,2)-2*size(os_10.est_epsilon,2)-1;
    gamma_gravity_2010 = os_10.ests2((end-start):end);

    % Baseline (taking the first stage results as given)
    [~,os_10] =  GMM_wrapper_gravity( ps_10, kappa_tau, kappa_epsilon,'base',0,gamma_gravity_2010,type );
    os_10.title = '2010 Data';     
    for t = 1:1:size(types,2)
        graph_overlay([ o os_10],[ d ps_10],types{t},'../../Output/F5_y2010',  "", 0)
    end
    


%% 2014
    fix = 2;
    knots = 3;
    M2014 = csvread(strcat(data_dir,'2014_Tetiuw.csv'));
    type = 'all';
    ps_14 = make_data_combo(M2014,knots,'cubic','','',0,type);    ps_14.sigma = sigma;       ps_14.x = x;
    [~,os_14] =  GMM_wrapper_gravity( ps_14, kappa_tau, kappa_epsilon,'base',fix,gamma_guess,type );
    start = size(os_14.ests2,2)-2*size(os_14.est_epsilon,2)-1;
    gamma_gravity_2014 = os_14.ests2((end-start):end);

    % Baseline (taking the first stage results as given)
    [~,os_14] =  GMM_wrapper_gravity( ps_14, kappa_tau, kappa_epsilon,'base',0,gamma_gravity_2014,type );
    os_14.title = '2014 Data';     
    for t = 1:1:size(types,2)
        graph_overlay([ o os_14],[ d ps_14],types{t},'../../Output/F5_y2014',  "", 0)
    end


%% Base 2012
    fix = 2;
    knots = 3;
    B2012 = csvread(strcat(data_dir,'2012_Tuw.csv'));
    type = 'all';
    ps_12B = make_data_combo(B2012,knots,'cubic','','',0,type);    ps_12B.sigma = sigma;       ps_12B.x = x;
    [~,os_12B] =  GMM_wrapper_gravity( ps_12B, kappa_tau, kappa_epsilon,'base',fix,gamma_guess,type );
    start = size(os_12B.ests2,2)-2*size(os_12B.est_epsilon,2)-1;
    gamma_gravity_2012T = os_12B.ests2((end-start):end);

    % Baseline (taking the first stage results as given)
    [~,os_12B] =  GMM_wrapper_gravity( ps_12B, kappa_tau, kappa_epsilon,'base',0,gamma_gravity_2012T,type );
    os_12B.title = '2012 Raw Tariffs';     
    for t = 1:1:size(types,2)
        graph_overlay([ o os_12B],[ d ps_12B],types{t},'../../Output/F5_T2012_v_TETI',  "", 0)
    end





    


%% HS
    basefile_Sf = '2012_Tuw_h0.csv'; name_Sf = '2012_T_h0';
    M_Sf = csvread(strcat(data_dir,basefile_Sf));
    knot_Sf = 3;
    fix = 2; % Fix = 2: means allow gamma to be estimated internally
    base_Sf = make_data_combo_hs_multi(M_Sf,knot_Sf,'cubic','','',0,'all' );     base_Sf.sigma = sigma;       base_Sf.x = x;
    gamma_guessHS = [0.5533   -0.0410   -0.0557   -0.4006   -0.3545    0.2085   -0.0830   -0.0148   -0.0016    0.06080];
    [~,obase_Sf] =  GMM_wrapper_gravity( base_Sf, kappa_tau, kappa_epsilon,'base',fix,gamma_guessHS,'all' );
    start = size(obase_Sf.ests2,2)-2*size(obase_Sf.est_epsilon,2)-1;
    gamma_gravity_HS = obase_Sf.ests2((end-start):end);
    Figure_Sequencer(base_Sf,obase_Sf,'F4_HS_all','')

    obase_Sf.title = "Sectoral-Level";
    base_Sf.kG = {};
    base_Sf = rmfield(base_Sf,'h');
    graph_overlay([ o obase_Sf],[ d base_Sf],'theta','../../Output/F4_HS_all',  "", 0)
    graph_overlay([ o obase_Sf],[ d base_Sf],'epsilon','../../Output/F4_HS_all',  "", 0)
    graph_overlay([ o obase_Sf],[ d base_Sf],'rho','../../Output/F4_HS_all',  "", 0)

% Wealth-Origin - HS Level
% Wealth Origin-Destination
   basefile_Sf = '2012_Tuw_h0.csv'; name_Sf = '2012_T_h0';
    
    M_Sf = csvread(strcat(data_dir,basefile_Sf));
    fix = 0;
    knot_S = 3;

    d2q_Wo_S = make_data_combo_hs_multi(M_Sf,knot_S,'cubic','wealth_oI','',0,'all' ); 
    d2q_Wo_S.sigma = sigma;  d2q_Wo_S.x = x;
    [~,o2q_Wo_S] =  GMM_wrapper_gravity( d2q_Wo_S, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity_HS,'all');
    Figure_Sequencer(d2q_Wo_S,o2q_Wo_S,'F4b_HS_origin','_split')


    d2q_Wd_S = make_data_combo_hs_multi(M_Sf,knot_S,'cubic','wealth_dI','',0,'all') ; 
    d2q_Wd_S.sigma = sigma;  d2q_Wd_S.x = x;
    [~,o2q_Wd_S] =  GMM_wrapper_gravity( d2q_Wd_S, kappa_tau, kappa_epsilon,'test',fix ,gamma_gravity_HS,'all');
    Figure_Sequencer(d2q_Wd_S,o2q_Wd_S,'F4b_HS_destination','_split')



%% HS
    hs = 3 ;  hs_routine_gravity; %  Foodstuffs
    hs = 4 ;  hs_routine_gravity; %  Mineral Products
    hs = 5 ;  hs_routine_gravity; %  Chemicals & Allied Industries
    hs = 6 ;  hs_routine_gravity; %  Plastics / Rubbers
    hs = 8 ;  hs_routine_gravity; %  Wood & Wood Products
    hs = 9 ;  hs_routine_gravity; %  Textiles
    hs = 12;  hs_routine_gravity; %  Metals
    hs = 13;  hs_routine_gravity; %  Machinery / Electrical
    hs = 14;  hs_routine_gravity; %  Transportation
    