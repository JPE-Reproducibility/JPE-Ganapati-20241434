function [ OUT,o  ] = GMM_wrapper_gravity( d, kappa_tau, kappa_epsilon,filename,fix,gamma,type )
%% Full GMM Wrapper
disp("Full GMM Wrapper")

% Fix = 0 ; base version, gamma COMPUTED INTERNALLY in an intermediate stage
% Fix = 1 ; old depreciated code could be recoded. TBA
% Fix = 2 ; fixed gamma - initial guess used in estimation


% SET UP DIMENTIONS
    R1 = 1:size(d.X_K,2);
    R2 = size(d.X_K,2)+1:2*size(d.X_K,2);
    R3 = max(R2)+1;

    switch type
        case 'all'
            gravity = [d.ldist_raw d.fta_wto d.comcur d.comlang d.colony ];
            if fix == 2
                gamma = [.3  0 0 0 0 .3 0 0 0 0];
            end
        case 'exog'
            gravity = [d.ldist_raw d.comcur d.comlang d.colony ];
            if fix == 2
                gamma = [.3  0 0 0  .3 0 0  0];
            end
        case 'dist'
            gravity = [d.ldist_raw ];
            if fix == 2
                gamma = [.3   .3 ];
            end
        case 'partial'
    end

%% First Stage Estimates
disp("First Stage Estimates")
%  Use identity weighting matrix
%  On-advice from CC - use identity weights only
    W_init = 1;
    tic
    if fix == 0
        init    = [R1*0+1 R2*0+1];
        C       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,0, W_init,fix,gamma,gravity );
        [ests]  = fminunc(@(x) C(x)   ,init,optimset('TolFun',1e-10,'TolX',1e-8,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
    elseif fix == 1

    elseif fix == 2
        init    = [R1*0+1 R2*0+1 gamma];
        C       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,0, W_init,fix,x(R3:end),gravity );
        [ests]  = fminunc(@(x) C(x)   ,init,optimset('TolFun',1e-10,'TolX',1e-8,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
        gamma = ests(R3:end);
    end
    toc

%% Create Optimal Second Stage Weighting Matrix
    disp("Create Optimal Second Stage Weighting Matrix")
    tic
    U  = GMM_gravity( d, ests(R1),ests(R2),kappa_tau,kappa_epsilon,1, W_init,fix,gamma,gravity );
    DZ_K = size(d.ZA_K,2);
    LambdaA = (d.ZA_K'.*repmat(U.G_A,1,DZ_K)')*(repmat(U.G_A,1,DZ_K).*d.ZA_K);
    LambdaB = (d.ZB_K'.*repmat(U.G_B,1,DZ_K)')*(repmat(U.G_B,1,DZ_K).*d.ZB_K);
    toc

% % Loop Manually      
%     LambdaA = zeros(size(d.ZA_K,2));
%     LambdaB = zeros(size(d.ZB_K,2));
%     tic
%     for i=1:size(U.G_A,1)
%         LambdaA = LambdaA +d.ZA_K(i,:)'*U.G_A(i)*U.G_A(i)'*d.ZA_K(i,:);
%         LambdaB = LambdaB +d.ZB_K(i,:)'*U.G_B(i)*U.G_B(i)'*d.ZB_K(i,:);
%     end
%     toc
%     sum(sum(abs(LambdaB-LambdaB_a)))

    tic
    Lambda  = blkdiag(LambdaA, LambdaB);
    Lambda  = Lambda/d.N;
    disp("Create Optimal Second Stage Weighting Matrix - Invert")
    W       = pinv(Lambda);
    toc

%% Second Stage
    disp("Second Stage")
    tic
    if fix == 0
        C       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,0, W_init,fix,gamma,gravity );
        [ests2]  = fminunc(@(x) C(x)   ,ests,optimset('TolFun',1e-10,'TolX',1e-12,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
        ests2 = [ests2 gamma];
    elseif fix == 1

    elseif fix == 2
        C       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,0, W_init,fix,x(R3:end),gravity );
        [ests2]  = fminunc(@(x) C(x)   ,ests,optimset('TolFun',1e-10,'TolX',1e-12,'MaxFunEvals',50000,'Display','iter','MaxIter',100000,'GradObj','off'));
    end
    toc

%% Standard Errors Use heteroskedastic SE
% Use perturbation method to get underlying matrices
% Do for just the spline, not the unreported nuisance paramters

    disp("Standard Errors Computing")
    if fix == 2
    estlin = [ests2(R1) ests2(R2) gamma];
    else
    estlin = [ests2(R1) ests2(R2)];
    end

    G_barA  = zeros(size(d.ZA_K,2),size(estlin,2));
    G_barB  = zeros(size(d.ZB_K,2),size(estlin,2));
    G       = @(x) GMM_gravity( d, x(R1),x(R2),kappa_tau,kappa_epsilon,2, W,fix,x(R3:end),gravity);
    u_diff = zeros(size(G(ests2),1),size(estlin,2));
    disp("Standard Errors Computing - Loopin 1")
    tic
    for dim = 1:size(estlin,2) 
        TOL             = .000001;
        pertub          = ests2;
        pertub(dim)     = pertub(dim) + TOL;
        u_diff(:,dim)   = (G(ests2)-G(pertub))/TOL;
    end
    toc
    tic
    for i=1:(size(d.ZA_K,1))
        G_barA          = G_barA +d.ZA_K(i,:)'*u_diff(i,:);
        G_barB          = G_barB +d.ZB_K(i,:)'*u_diff(i+d.N,:);
    end
    G_bar   = [G_barA ; G_barB]/d.N;
    toc

    disp("Standard Errors Computing - Loopin 2")
    U       = GMM_gravity( d, ests2(R1),ests2(R2),kappa_tau,kappa_epsilon,1, W,fix ,ests2(R3:end),gravity);
    tic
    LambdaA = (d.ZA_K'.*repmat(U.G_A,1,DZ_K)')*(repmat(U.G_A,1,DZ_K).*d.ZA_K);
    LambdaB = (d.ZB_K'.*repmat(U.G_B,1,DZ_K)')*(repmat(U.G_B,1,DZ_K).*d.ZB_K);
    LambdaAB = (d.ZA_K'.*repmat(U.G_A,1,DZ_K)')*(repmat(U.G_B,1,DZ_K).*d.ZB_K);
    toc

% % Loop Manually      
%     LambdaA = zeros(size(d.ZA_K,2));
%     LambdaB = zeros(size(d.ZB_K,2));
%     LambdaAB= zeros(size(LambdaA,1),size(LambdaB,1));
%     tic
%     d.ZA_K = sparse(d.ZA_K);
%     U.G_A = sparse(U.G_A);
%     d.ZB_K = sparse(d.ZB_K);
%     U.G_B = sparse(U.G_B);
%     for i=1:size(U.G_A,1)
%         LambdaA = LambdaA +d.ZA_K(i,:)'*U.G_A(i)*U.G_A(i)'*d.ZA_K(i,:);
%         LambdaB = LambdaB +d.ZB_K(i,:)'*U.G_B(i)*U.G_B(i)'*d.ZB_K(i,:);
%         LambdaAB = LambdaAB +d.ZA_K(i,:)'*U.G_A(i)*U.G_B(i)'*d.ZB_K(i,:);
%     end
%     toc
%     sum(sum(abs(LambdaAB-LambdaAB_a)))

    disp("Standard Errors Computing - Inverting")
    tic
    Lambda = [LambdaA LambdaAB; LambdaAB' LambdaB];
    Lambda  = Lambda/d.N;
    W       = pinv(Lambda);
    toc
    tic
    V       = ((G_bar'*W*G_bar)^-1*G_bar'*W*Lambda*W*G_bar*(G_bar'*W*G_bar)^-1)/d.N;
    SE      = diag(V).^.5;
    OUT     = [ests2(1:size(estlin,2)); SE'];
    toc


%% Standard Errors full functional form
    disp("Code From CC 2017 QE")
        try
            [o.standard_errors1] = CC_se(R1,ests2,d,U.G_A,G_bar,W,Lambda,V);
            [o.standard_errors2] = CC_se(R2,ests2,d,U.G_B,G_bar,W,Lambda,V);
%             [o.standard_errors_theta_e] = CC_se_thetae(R1,R2,ests2,d,U.G_A,U.G_B,G_bar,W,Lambda,V);
            [o.standard_errors_theta_e] = CC_se_thetae(R1,R2,ests2,d,U.G_A,U.G_B,G_bar,W,Lambda,V);

        catch
            disp("CC Standard Errors Failed")
        end

%% Output Data
    disp("Output Data")
    o.G_bar     = G_bar;
    o.Lambda    = Lambda;
    o.W         = W;
    o.SE        = SE;
    o.V         = V;
    o.ests2     = ests2;
    o.est_epsilon=ests2(R1);
    o.est_rho    =ests2(R2);
    o.filename  = filename;

end