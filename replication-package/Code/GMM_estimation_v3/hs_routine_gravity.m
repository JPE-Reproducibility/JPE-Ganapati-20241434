clc

% Use Gamma from ALL DATA
fix = 0;
kappa_tau = 1;
kappa_f   = 0;
kappa_epsilon = 1/((sigma-1)*kappa_tau+kappa_f);

    switch hs
      case 1 
        title = "Animal & Animal Products";
      case 2 
        title = "Vegetable Products";
      case 3 
        title = "Foodstuffs";
      case 4 
        title = "Mineral Products";
      case 5 
        title = "Chemicals & Allied Industries";
      case 6 
        title = "Plastics / Rubbers";
      case 7 
        title = "Raw Hides, Skins, Leather, & Furs";
      case 8 
        title = "Wood & Wood Products";
      case 9 
        title = "Textiles and Apparel";
      case 10
        title = "Footwear / Headgear";
      case 11
        title = "Stone / Glass";
      case 12
        title = "Metals";
      case 13
        title = "Machinery / Electrical";
      case 14
        title = "Transportation";
      case 15
        title = "Misc Equipment";
      case 16
        title = "Miscellaneous";
      case 99
        title = "composite";
      otherwise 
        title = "misc";
    end

%% Baseline
    M_hs = csvread(strcat(data_dir,['2012_Tuw_h' int2str(hs) '.csv']));
    knot_hs = 33;
    fix = 0;  %% THIS LOCKS IT INTO PLACE from that for linear one
    d_hs = make_data_combo(M_hs,knot_hs,'cubic','','',0,'dist');     d_hs.sigma = sigma;       d_hs.x = x;
    [~,o_hs] =  GMM_wrapper_gravity( d_hs, kappa_tau, kappa_epsilon,['base'],fix,gamma_gravity_HS,'all');
    o_hs.title = 'Semi-parametric';     
    p_overlay_elasticity( d_hs,o_hs,'theta'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat(title))
    p_overlay_elasticity( d_hs,o_hs,'epsilon'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat('Extensive: ',title))
    p_overlay_elasticity( d_hs,o_hs,'rho'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat(title))
    p_overlay_elasticity( d_hs,o_hs,'extensive'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat(title))

    
%% Wealth Origin-Destination
    fix = 0;
    d_split_hs = make_data_combo(M_hs,knot_hs,'cubic','wealth_oI','',0,'dist'); 
    d_split_hs.sigma = sigma;  d_split_hs.x = x;
    [~,o_split_hs] =  GMM_wrapper_gravity( d_split_hs, kappa_tau, kappa_epsilon,['test'],fix ,gamma_gravity_HS,'all');
    p_overlay_elasticity( d_split_hs,o_split_hs,'theta_split'       ,['../../Output/F6_HS_' int2str(hs)] ,strcat('\theta: ',title))

