function [ OUT,o  ] = GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,~,type )
%% Full GMM Wrapper
disp("Full GMM Wrapper")

% SET UP DIMENTIONS
    R1 = 1:size(d.X_K,2);
    R2 = size(d.X_K,2)+1:2*size(d.X_K,2);

%% First Stage Estimates
disp("First Stage Estimates")
%  Use identity weighting matrix
%  On-advice from CC - use identity weights only
    W_init = 1;
    tic
    init    = [R1*0 R2*0];
    % histogram(d.R_nE_ij)

    if strcmp(type,"LogNormal")
        init= .76;
        C       = @(x) GMM_gravity_log_normal( d, 0,x,kappa_tau,kappa_epsilon,0, W_init );
        [ests]  = fminunc(@(x) C(x)   ,init,optimset('TolFun',1e-10,'TolX',1e-8,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
        ests = [0 ests];
    else
        init(1) = -.8;
        init(size(d.X_K,2)+1)= -.8;

        C       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,0, W_init );
        [ests]  = fminunc(@(x) C(x)   ,init,optimset('TolFun',1e-10,'TolX',1e-8,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
    end

    toc

%% Create Optimal Second Stage Weighting Matrix
    disp("Create Optimal Second Stage Weighting Matrix")
    U  = GMM_gravity( d, ests(R1),ests(R2),kappa_tau,kappa_epsilon,1, W_init );
    DZ_K = size(d.ZA_K,2);
    LambdaA = (d.ZA_K'.*repmat(U.G_A,1,DZ_K)')*(repmat(U.G_A,1,DZ_K).*d.ZA_K);
    LambdaB = (d.ZB_K'.*repmat(U.G_B,1,DZ_K)')*(repmat(U.G_B,1,DZ_K).*d.ZB_K);
    toc

    Lambda  = blkdiag(LambdaA, LambdaB);
    Lambda  = Lambda/d.N;
    disp("Create Optimal Second Stage Weighting Matrix - Invert")
    W       = pinv(Lambda);
    toc

%% Second Stage
    disp("Second Stage")
    if strcmp(type,"LogNormal")
        C       = @(x) GMM_gravity_log_normal( d, 0,x,kappa_tau,kappa_epsilon,0, W_init );
        [ests2]  = fminunc(@(x) C(x)   ,ests(2),optimset('TolFun',1e-10,'TolX',1e-12,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
        ests2 = [0 ests2];
    else
        C       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,0, W_init );
        [ests2]  = fminunc(@(x) C(x)   ,ests,optimset('TolFun',1e-10,'TolX',1e-12,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
    end
    toc

%% Standard Errors Use heteroskedastic SE
% Use perturbation method to get underlying matrices
% Do for just the spline, not the unreported nuisance paramters

    disp("Standard Errors Computing")
    estlin = [ests2(R1) ests2(R2)];

    G_barA  = zeros(size(d.ZA_K,2),size(estlin,2));
    G_barB  = zeros(size(d.ZB_K,2),size(estlin,2));
    if strcmp(type,"LogNormal")
        G       = @(x) GMM_gravity_log_normal( d, x(R1),x(R2),kappa_tau,kappa_epsilon,2, W);
    else
        G       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,2, W);
    end
    u_diff = zeros(size(G(ests2),1),size(estlin,2));
    disp("Standard Errors Computing - Loopin 1")
    for dim = 1:size(estlin,2) 
        TOL             = .000001;
        pertub          = ests2;
        pertub(dim)     = pertub(dim) + TOL;
        u_diff(:,dim)   = (G(ests2)-G(pertub))/TOL;
    end

    for i=1:(size(d.ZA_K,1))
        G_barA          = G_barA +d.ZA_K(i,:)'*u_diff(i,:);
        G_barB          = G_barB +d.ZB_K(i,:)'*u_diff(i+d.N,:);
    end
    G_bar   = [G_barA ; G_barB]/d.N;
    toc

    disp("Standard Errors Computing - Loopin 2")
    if strcmp(type,"LogNormal")
        U       = GMM_gravity_log_normal( d, ests2(R1),ests2(R2),kappa_tau,kappa_epsilon,1, W);
    else
        U       = GMM_gravity( d, ests2(R1),ests2(R2),kappa_tau,kappa_epsilon,1, W);
    end
    LambdaA = (d.ZA_K'.*repmat(U.G_A,1,DZ_K)')*(repmat(U.G_A,1,DZ_K).*d.ZA_K);
    LambdaB = (d.ZB_K'.*repmat(U.G_B,1,DZ_K)')*(repmat(U.G_B,1,DZ_K).*d.ZB_K);
    LambdaAB = (d.ZA_K'.*repmat(U.G_A,1,DZ_K)')*(repmat(U.G_B,1,DZ_K).*d.ZB_K);
    toc

    disp("Standard Errors Computing - Inverting")
    Lambda = [LambdaA LambdaAB; LambdaAB' LambdaB];
    Lambda  = Lambda/d.N;
    W       = pinv(Lambda);
    V       = ((G_bar'*W*G_bar)^-1*G_bar'*W*Lambda*W*G_bar*(G_bar'*W*G_bar)^-1)/d.N;
    SE      = diag(V).^.5;
    OUT     = [ests2(1:size(estlin,2)); SE'];
    toc

%% Output Data
    disp("Output Data")
    %     o.G_bar     = G_bar;
    %     o.Lambda    = Lambda;
    %     o.W         = W;
    o.SE        = SE;
    o.V         = V;
    o.ests2     = ests2;
    o.knots     = d.kG;
    if strcmp(type,"LogNormal")
        o.est_sigma = ests2(R2);
    else
        o.est_epsilon=ests2(R1);
        o.est_rho    =ests2(R2);
    end
end