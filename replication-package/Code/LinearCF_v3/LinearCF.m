function [Real_Wage,dN_linear,labels,trade_share,out] = LinearCF( segments,sigma,scenario,spline)

file_base = '2012B.csv';
addpath([  '../GMM_estimation_v3'])
data_dir = [ '../../Data/Int/WIOD_sampleB/'];

testing  = 0;
balanced = 0;

% Code for testing
if testing == 1
    scenario    = 0;
    sigma       = 4;
    balanced    = 1;
end


%% Load in data
% notation: first  index is ALWAYS i
%           second index is always j
%           third  index is sequence (when needed)
disp('Loading Data');

if testing == 1
    RandStream.setGlobalStream(RandStream('mcg16807','Seed',14500))
    dim         = 2;
    n_ij        = rand(dim);
    n_ij        = diag(.9)+(.25*(n_ij-diag(diag(n_ij))));
    xbar_ij     = 1/2*(rand(dim)+n_ij)+1/2;
    N_i         = 100+0*rand(dim,1);

    % make symmetrical.
    xbar_ij     = (xbar_ij + xbar_ij')/2;
    n_ij        = (n_ij + n_ij')/2;
    G_ij        = randi(4,dim,dim);

    labels = 1;
    
else
    % n_ij        = csvread([ ROOT '/Data/Int/WIOD_sampleB/n_ij_'     file_base ]);
    % X_ij        = csvread([ ROOT '/Data/Int/WIOD_sampleB/XX_ij_'    file_base ]);
    % N_ii        = csvread([ ROOT '/Data/Int/WIOD_sampleB/N_ii_'     file_base ]);
    % G_ij        = csvread([ ROOT '/Data/Int/WIOD_sampleB/G_ij_'     file_base ]);
    % labels      = importdata([ ROOT '/Data/Int/WIOD_sampleB/Names_' file_base ]);
    n_ij        = csvread([ data_dir 'n_ij_'     file_base ]);
    X_ij        = csvread([ data_dir 'XX_ij_'    file_base ]);
    N_ii        = csvread([ data_dir 'N_ii_'     file_base ]);
    G_ij        = csvread([ data_dir 'G_ij_'     file_base ]);
    labels      = importdata([ data_dir 'Names_' file_base ]);

    dim         = size(n_ij,1);
    N_i         = N_ii(:,1)./diag(n_ij);
    N_ij        = n_ij.*repmat(N_i,1,dim);
    xbar_ij     = X_ij./(N_ij);
end


%% Create Share data (lowercase x and y)
X_ij = xbar_ij.*n_ij.*repmat(N_i,1,dim);
out.X_ij_init = X_ij;

if testing == 1
    X_ij = [.15 .35; .35 .15];
end

% iota - equation 9 dlnE=i_j.*dlnw_i E = w*L+T
% first is vector; second is matrix on i dimention
i_vec = sum(X_ij,2)./ sum(X_ij,1)';
i_i = repmat(sum(X_ij,2)./ sum(X_ij,1)',1,dim);

if balanced == 1
    % % This is to set trade to be balanced
    factor = (sum(X_ij,2) - sum(X_ij,1)')/dim;
    sum(X_ij - factor,2)./ sum(X_ij - factor,1)';
    X_ij = X_ij - factor;
    
    % Normalize world income to 1
    X_ij = X_ij./sum(sum(X_ij));
    i_vec = ones(size(i_vec)); 
    i_i = ones(size(i_i));
end

% Create Matrices
x_ij = X_ij ./ repmat(sum(X_ij,1),dim,1);
y_ij = X_ij ./ repmat(sum(X_ij,2),1,dim);
trade_share = sum(y_ij,2)-diag(y_ij);


%% Code in Elasticities
if testing == 1
    elast_r         = @(x) -.7.*1+x*0;
    elast_e         = @(x) -.7.*1+x*0;
    theta           = @(n) (1-sigma).*(1+elast_r(n)  - elast_e(n) )./elast_e(n)  ;
    elast_rho       = elast_r;
    elast_epsilon   = elast_e;

    Gelast_epsilon   = @(x,G) elast_epsilon(x);
    Gelast_rho       = @(x,G) elast_rho(x);
    Gtheta           = @(n,G) (1-sigma)*(1+elast_rho(n) - elast_epsilon(n))./elast_epsilon(n);

elseif testing ~= 1

    kappa_tau = 1;
    kappa_f   = 0;
    sigma = 3.2;
    kappa_epsilon = 1/((sigma-1)*kappa_tau+kappa_f);
%    (from Teti Tariffs) gamma_gravity = [ 0.3520   -0.0214   -0.0189   -0.2324   -0.1502    0.2812   -0.0834   -0.0015   -0.0576   -0.0471];
%    (from Base Tariffs) gamma_gravity = [ 0.3920   -0.0408   -0.0131   -0.2596   -0.1635    0.3402   -0.1012    0.0039   -0.0985   -0.0671];
    load('../GMM_estimation_v3/gamma_gravity.mat');

    fix = 0;    
    space           = .0005;
    unif            = 0:space:(1-space);
    range           = 1:1:1000;
    x               = log(unif(range));

    basefile = '2012_Tetiuw.csv';
    M = csvread(strcat(data_dir,basefile));

    if spline == 0
        if isfile('ESTIMATES_LINEAR.mat')
            disp('Loading estimates')
            load('ESTIMATES_LINEAR.mat');
        else
            d = make_data_combo(M,1,'cubic','','',0,'all');
            d.sigma = sigma;  d.x = x;
            [~,o] =  GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,'CF_SPLINE',fix ,gamma_gravity,'all');

            save('ESTIMATES_LINEAR','d','o','G_ij')
        end
            G_ij = zeros(size(X_ij));
    elseif spline == 2
        if isfile('ESTIMATES_SPLINE2.mat')
            disp('Loading estimates')
            load('ESTIMATES_SPLINE2.mat');
        else

            d = make_data_combo(M,3,'cubic','','',0,'all');
            d.sigma = sigma;  d.x = x;
            [~,o] =  GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,'CF_SPLINE',fix ,gamma_gravity,'all');
            save('ESTIMATES_SPLINE2','d','o','G_ij')
        end
        G_ijr        = csvread([  '../../Data/Int/WIOD_sampleB/G_ij_'     file_base ]);
        G_ij = zeros(size(G_ijr));
      
    elseif spline == 4
        if isfile('ESTIMATES_SPLINE4.mat')
            disp('Loading estimates')
            load('ESTIMATES_SPLINE4.mat');
        else

            d = make_data_combo(M,3,'cubic','wealth_oI','',0,'all');
            d.sigma = sigma;  d.x = x;
            [~,o] =  GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,'CF_SPLINE',fix ,gamma_gravity,'all');
            %p_overlay_elasticity( d,o,'rho_split'       ,'test' );
            save('ESTIMATES_SPLINE4','d','o','G_ij')
        end
        G_ijr        = csvread([  '../../Data/Int/WIOD_sampleB/G_ij_'     file_base ]);
        G_ij = zeros(size(G_ijr));
        G_ij(G_ijr == 0) = 1;
        G_ij(G_ijr == 1) = 1;
        G_ij(G_ijr == 2) = 0;
        G_ij(G_ijr == 3) = 0;
    elseif spline == 6
        if isfile('ESTIMATES_SPLINE6.mat')
            disp('Loading estimates')
            load('ESTIMATES_SPLINE6.mat');
        else

            d = make_data_combo(M,3,'cubic','wealth_odI','',0,'all');
            d.sigma = sigma;  d.x = x;
            [~,o] =  GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,'CF_SPLINE',fix ,gamma_gravity,'all');
            %p_overlay_elasticity( d,o,'rho_split'       ,'test' );
            save('ESTIMATES_SPLINE6','d','o','G_ij')
        end
        G_ijr        = csvread([  '../../Data/Int/WIOD_sampleB/G_ij_'     file_base ]);
        G_ij = G_ijr;
    end

    GDerivitive      = @(n,G) evknots_make_derivG(d.kG,log(n),G) ;
    Gelast_epsilon   = @(x,G) reshape( GDerivitive((x(:)),G(:))*o.est_epsilon',dim,dim);
    Gelast_rho       = @(x,G) reshape( GDerivitive((x(:)),G(:))*o.est_rho',dim,dim);
    Gtheta           = @(n,G) (1-sigma)*(1+Gelast_rho(n,G) - Gelast_epsilon(n,G))./Gelast_epsilon(n,G);

end



%% Set Up Trade Segments
tic
if testing == 1
    index= 1;
    r_1           = ones(dim,dim);     r_2           = r_1*exp(0.2);
    r_1(logical(eye(size(r_1)))) = 1;  r_2(logical(eye(size(r_2)))) = 1;
    f_1           = ones(dim,dim);     f_2           = ones(dim,dim); 

end

y_ij_init = y_ij;
x_ij_init = x_ij;
X_ij_init = X_ij;
n_ij_init = n_ij;

F_1   = ones(dim,1);    
F_2   = ones(dim,1);
L_1   = ones(dim,1);    
L_2   = ones(dim,1);
T_1   = ones(dim,1);    
T_2   = ones(dim,1);

if scenario == 2

    r_1           = ones(dim,dim);
    r_2           = (r_1*.99).^(1/(1-sigma));

    f_1           = ones(dim,dim);
    f_2           = f_1*1;

    r_1(logical(eye(size(r_1)))) = 1;    r_2(logical(eye(size(r_2)))) = 1;
    f_1(logical(eye(size(f_1)))) = 1;    f_2(logical(eye(size(f_1)))) = 1;

    out.file = '1pc';
    out.title= '1% change in r_ij';
    out.shock =.01;

elseif scenario == 5

    r_1           = ones(dim,dim);
    G_ij_orig     = csvread([ data_dir 'G_ij_'      file_base ]);
    G4_ij_orig    = csvread([ data_dir 'G4_ij_'     file_base ]);
    GSP_ij_delta_ = csvread([ data_dir 'GSP_delta_' file_base ]);
    GSP_ij = GSP_ij_delta_>1;

 
    out.GSP_ij      = GSP_ij;
    out.G_ij_orig   = G_ij_orig;
    out.G4_ij_orig  = G4_ij_orig;
    
    % Only to Developing Countries
    r_2             = ones(dim,dim);
    r_2(GSP_ij==1)  = (.99).^(1/(1-sigma));
    f_1             = ones(dim,dim);
    f_2             = f_1*1;

    r_1(logical(eye(size(r_1)))) = 1;    r_2(logical(eye(size(r_2)))) = 1;
    f_1(logical(eye(size(f_1)))) = 1;    f_2(logical(eye(size(f_1)))) = 1;

    out.file    = 'GSP_like';
    out.title   = 'GSP Trade Facilitation';
    out.shock   = .1;

end


%% Setup interpolation
    if min(min(f_1-f_2)) == 0 && min(max(f_1-f_2)) == 0
    else
        % We need to play with the lower bounds on both this and the zero
        buffer = 1e-11;
        x = linspace(.000001,1,2000);

        Val_interp = [];
        for g = 1:max(size(d.kG))
            GratioA0 = @(nn) GratioA(nn,(g-1));
            GrhoA0 = @(nn) GrhoA(nn,(g-1));   
            Gepsilon_barA0 = @(nn) Gepsilon_barA(nn,(g-1)); 
            topf0 = arrayfun(@(lb,ub) integral(GratioA0,lb,ub) , zeros(1,2000)+buffer,x);
            botf0 = arrayfun(@(lb,ub) integral(GrhoA0,lb,ub)  , zeros(1,2000)+buffer,x)./Gepsilon_barA0(x);
            i0 = topf0./botf0;
            Val_interp = [Val_interp; i0];
        end

        range = 0:1:(max(size(d.kG))-1);

        if max(size(d.kG)) > 1
            Gpi_inter = griddedInterpolant(repmat(x,max(size(d.kG)),1)',repmat(range',1,2000)',Val_interp');
        else
            pi_inter = griddedInterpolant(repmat(x,max(size(d.kG)),1)',Val_interp');
            Gpi_inter = @(n,G) pi_inter(n);
        end

    end

% %% Initialize Segments
lr_segments      = Create_Segments(log(r_1),log(r_2),segments);
lf_segments      = Create_Segments(log(f_1),log(f_2),segments);
lF_segments      = Create_Segments(log(F_1),log(F_2),segments);
lL_segments      = Create_Segments(log(L_1),log(L_2),segments);
lT_segments      = Create_Segments(log(T_1),log(T_2),segments);

%% Initialize Loop
dlnw_sequence   = zeros(dim,segments);
dlnP_sequence   = zeros(dim,segments);
dlnN_sequence   = zeros(dim,segments);
dlnn_sequence   = zeros(dim,dim,segments);
dlnxbar_sequence= zeros(dim,dim,segments);

%% Do Loop
for count = 1:segments

%     disp(count)
    dln_r = (lr_segments(:,:,count+1))-(lr_segments(:,:,count));
    dln_f = (lf_segments(:,:,count+1))-(lf_segments(:,:,count));
    dln_F = (lF_segments(:,:,count+1))-(lF_segments(:,:,count));
    dln_T = (lT_segments(:,:,count+1))-(lT_segments(:,:,count));
    dln_L = (lL_segments(:,:,count+1))-(lL_segments(:,:,count));

    theta_ij = Gtheta(n_ij,G_ij);
    theta_ijS = theta_ij./(sigma-1);

    % Shifter Functions
    dln_r_p = sum(x_ij.* ...
        (theta_ij./(1-sigma).*dln_r ...
        - repmat(sum(y_ij.*theta_ij./(1-sigma).*dln_r,2),1,dim)) ...
        ,1)';

    dln_f_p = sum(x_ij.* ...
        ((1-theta_ijS).*dln_f ...
        - repmat(sum(y_ij.*(1-theta_ijS).*dln_f,2),1,dim)) ...
        ,1)';

    dln_L_p = (1-sum(x_ij.*theta_ijS,1)).*i_vec'.*dln_L' ...
        - sum(x_ij.*repmat(dln_L ...
        - sum(y_ij.*theta_ijS.*repmat((i_vec.*dln_L)',dim,1),2),1,dim) ...
        ,1);
    dln_L_p = dln_L_p';

    dln_T_p = (1-sum(x_ij.*theta_ijS,1)) .*(1-i_vec)'.*dln_T' ...
        + sum(x_ij.*repmat(sum(y_ij.*theta_ijS.*repmat(((1-i_vec).*dln_T)',dim,1),2),1,dim) ...
        ,1);
    dln_T_p = dln_T_p';


    % Wage Aux Functions
    v_w = zeros(dim);
    v_p = zeros(dim);

    % These are the most likley to contain an error
    theta_calc = Gtheta(n_ij,G_ij);
    for i=1:dim
        for j = 1:dim
            if i  == j
                temp = 0;
                for o = 1:dim
                    temp = temp + x_ij(o,i)*(theta_calc(o,i));
                    v_p(i,j) = v_p(i,j) + x_ij(o,i)*theta_calc(o,i);
                end
                v_w(i,j) = v_w(i,j) + (1-temp/(sigma-1))*i_vec(i);
            end

            for o = 1:dim
                v_w(i,j) = v_w(i,j)+x_ij(o,i)*y_ij(o,j)*(theta_calc(o,j))*i_vec(j)/(sigma-1);
                v_p(i,j) = v_p(i,j)-x_ij(o,i)*y_ij(o,j)*(theta_calc(o,j));
            end

            temp = 0;
            for dd = 1:dim
                temp = temp + y_ij(j,dd)*(theta_calc(j,dd));
            end
            v_w(i,j) = v_w(i,j)- x_ij(j,i)*(1-sigma/(sigma-1)*((theta_calc(j,i))-temp));
        end
    end

    % Price Aux Function
    m_w = diag(sigma)*diag(ones(dim,1)) -    y_ij.*   i_i';
    m_p = y_ij*(sigma-1);

    % Wage Equations
    dln_r_w = sum(y_ij.*dln_r,2);
        
    dln_T_w = sum(y_ij.*(1-i_i').*repmat(dln_T',dim,1),2);

    if min(min(dln_f)) == 0 && min(max(dln_f)) == 0
        dln_f_w = 0;
    else
        dln_f_w = sum(y_ij.*Gpi_inter(n_ij,G_ij).*dln_f,2);
    end
        
    if min(dln_F) == 0 && max(dln_F) == 0
        dln_F_w = 0;
        dln_L_w = 0;
    else
        % Pi
        pi_i = 1 - sum(y_ij.*Gpi_inter(n_ij,G_ij),2);
        dln_F_w = -pi_i.*dln_F;
        dln_L_w = -(1+pi_i).*dln_L + sum(y_ij.*i_i'.*repmat(dln_L',dim,1),2);
    end


    % Phi Equations
    dln_phi_p = dln_r_p - dln_f_p           + dln_L_p + dln_T_p;
    dln_phi_w = dln_r_w - dln_f_w + dln_F_w + dln_L_w + dln_T_w;

    % % Normalization Solver
    % This is slow and we really used use a more exact linear solver.
    % (as below), but that need some more debugging.
    FX = @(lw,lp) [dln_phi_p ; dln_phi_w] - [v_p*lp - v_w*lw ; -m_p*lp + m_w*lw];
    FXX = @(x) FX(x(1:dim,1),x((dim+1):2*dim,1));
    dlnw_solver2 = fsolve(@(x) FXX([0;x]),[zeros(dim-1,1);zeros(dim,1)],optimoptions('fsolve','Algorithm','levenberg-marquardt','FunctionTolerance',1e-15,'Display','off'));
    dlnw = [0; dlnw_solver2(1:dim-1)];

    %  solve dlnP 
    dlnP = dlnw_solver2(dim:end);

    dlnE = i_vec .* (dlnw + dln_L) + (1-i_vec).*dln_T;

    % solve dlnn, dlnx_bar EQ 48
    dlnn = (dln_f-dln_r ...
        + repmat(sigma*dlnw,1,dim) ...
        - repmat((sigma-1)*dlnP',dim,1) ...
        - repmat(dlnE',dim,1)) ...
    ./Gelast_epsilon(n_ij,G_ij)    ;

    %     solve dlnN EQ 53
    dlnN = dln_L - sum(y_ij.*dln_f,2)+sum(y_ij.*Gtheta(n_ij,G_ij).*Gelast_epsilon(n_ij,G_ij).*dlnn./(sigma-1),2);
    dlnx_bar = dln_f + repmat(dlnw,1,dim) + (Gelast_rho(n_ij,G_ij)-Gelast_epsilon(n_ij,G_ij)).*dlnn;

    %compute dlnX=dlnx_bar+dlnn+dlnN
    dlnX=dlnx_bar+dlnn+dlnN;

    %update lnX=lnX_prior+dlnX
    %update lnn=lnn_prior+dlnn
    X_ij =X_ij.*exp(dlnX);
    n_ij =n_ij.*exp(dlnn);

    %% Aggregate Up Welfare By Sequence
    TechShock  = 1/(sigma-1).*sum(x_ij.*(dln_r),1);
    ToT        = -sum(x_ij.*(repmat(dlnw,1,dim)-repmat(dlnw,1,dim)'),1);
    Extensive  = 1/(sigma-1).* sum(x_ij.*repmat((dlnN),1,dim),1);
    Selection  = 1/(sigma-1).* sum((Gelast_rho(n_ij_init,G_ij)+1).*x_ij.*((dlnn)),1);
    Error =   ((dlnw-dlnP) - TechShock' - ToT' - Extensive' -Selection')';
    Welfare = (dlnw-dlnP);

    TechShock_seq(:,count)  = TechShock';
    ToT_seq(:,count)        = ToT';
    Extensive_seq(:,count)  = Extensive';
    Selection_seq(:,count)  = Selection';
    Error_seq(:,count)      = Error';
    Welfare_seq(:,count)    = Welfare';

    segment_diag = [ dlnw-dlnP TechShock' ToT' Extensive' Selection' Error'];
%     segment_diag2 = mean(segment_diag)

    %    save sequence dlnw, dlnP, dlnN
    dlnw_sequence(:,count) = dlnw;
    dlnP_sequence(:,count) = dlnP;
    dlnN_sequence(:,count) = dlnN;
    dlnn_sequence(:,:,count) = dlnn;
    dlnxbar_sequence(:,:,count) = dlnx_bar;

    %compute x,y,i; epsilon,rho;theta
    x_ij = X_ij ./ repmat(sum(X_ij,1),dim,1);
    y_ij = X_ij ./ repmat(sum(X_ij,2),1,dim);
    i_i = repmat(sum(X_ij,2)./ sum(X_ij,1)',1,dim);
    i_vec = sum(X_ij,2)./ sum(X_ij,1)';

    if count == 1 && segments == 2
        %save('test','X_ij','n_ij','x_ij','y_ij','i_i','i_vec')
    end

    if sum(isnan(X_ij)) >0
        disp('ERROR')
        stop;
        break;
    end

end
toc

% Save Reults
dN_linear   = exp(sum(dlnN_sequence,2));
dw_linear   = exp(sum(dlnw_sequence,2));
dP_linear   = exp(sum(dlnP_sequence,2));
dn_linear   = exp(sum(dlnn_sequence,3));
dxbar_linear= exp(sum(dlnxbar_sequence,3));

Real_Wage   = dw_linear./dP_linear;
out.dN_linear = dN_linear;
out.dw_linear = dw_linear;
out.dP_linear = dP_linear;
out.dn_linear = dn_linear;
out.dxbar_linear = dxbar_linear;
out.Real_Wage = Real_Wage;
out.x_ij_hat = x_ij./x_ij_init;
out.f_ij_hat = f_2./f_1;
out.r_ij_hat = r_2./r_1;
out.F_i_hat  = F_2./F_1;
out.G_ij = G_ij;
out.i_vec = i_vec;
out.y_ij_init = y_ij_init;
out.x_ij_init = x_ij_init;
out.n_ij_init = n_ij_init;
out.X_ij = X_ij;
out.X_ij_init = X_ij_init;


%% Holds in All Decompositions

LRW = log(Real_Wage');
A_tech      = 1/(sigma-1).*sum(x_ij_init.*(log(r_2)-log(r_1)),1);
B_tot       = -sum(x_ij_init.*log(repmat(dw_linear,1,dim)./repmat(dw_linear,1,dim)'),1);
C_demand    = -1/(sigma-1).*sum(x_ij_init.*(log(x_ij)-log(x_ij_init)),1);
D_Extensive = 1/(sigma-1).* sum(x_ij_init.*repmat(log(dN_linear),1,dim),1);
E_Selection = 1/(sigma-1).* sum((Gelast_rho(n_ij_init,G_ij)+1).*x_ij_init.*(log(dn_linear)),1);
Error       =   (LRW' - A_tech' - B_tot' - C_demand' -D_Extensive' -E_Selection')';
out.Decomp_5 = [ LRW' A_tech' B_tot' C_demand' D_Extensive' E_Selection' Error'];


%% Test Code
Test_Code

end
