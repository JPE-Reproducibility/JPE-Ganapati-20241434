function [standard_errors] = CC_se_thetae(i_eps,i_rho,ests2,d,uhat_eps,uhat_rho,G_bar,W,Lambda,V)
%% Standard Errors full functional form
% Code From CC 2017 QE
        
%         [~,~,Deriv]=evknots_make(d.k,d.x);
%         uhat  = U.G_A;
        percentile_se = norminv([0.025 0.975]);

        [Psixx,~,Deriv] =evknots_make(d.k,d.x);
        Nboot = 1000;
        B = d.ZA_K;
        [n,~] = size(B);
        alpha = [ 0.9500 ];

        SGl   = inv(G_bar'*W*G_bar)*G_bar'*W;
        Q     = SGl*Lambda*SGl';

        nxx   = length(Psixx);
        Sxx_eps   = zeros(nxx,1);

        rbot = 1;
        rtop = size(SGl,2)/2;
        SGl_eps   = SGl(i_eps,rbot:rtop);
        Q_eps     = Q(i_eps,i_eps);

        %% For Derivative
        for i = 1:nxx
            Sxx_eps(i) = sqrt(6*Deriv(i,:)*Q_eps*Deriv(i,:)');
        end

        Deriv_eps = Deriv*ests2(i_eps)';
        invDeriv_eps = -1./(Deriv*ests2(i_eps)');

        %% GMM
        SE_EQ_TCC = zeros(size(Deriv,1),1);
        for i=1:size(Deriv,1)

            grad(i) = (1/Deriv_eps(i).^2);
            SE_EQ_TCC(i) = (grad(i)'*(Deriv(i,:)*V(i_eps ,i_eps )*Deriv(i,:)')*grad(i))^.5;

        end

        %% CC
        zb_eps   = zeros(Nboot,1);
        for i = 1:Nboot
            % Mammen two-point
            w     = 0.5 - 0.5*sqrt(5)*(2*(rand(n,1) < (sqrt(5)+1)/(2*sqrt(5)))-1);
            R_eps     = SGl_eps*(B'*(uhat_eps.*w))/sqrt(n);
            Rb_eps    = (Deriv*R_eps)./Sxx_eps;
            zb_eps(i) = max(abs(Rb_eps));
        end
        zalpha_eps = quantile(zb_eps,alpha);
        v_eps     = Sxx_eps/sqrt(n)*zalpha_eps.*(grad');

        Deriv_LB_CC    = invDeriv_eps - v_eps;
        Deriv_UB_CC    = invDeriv_eps + v_eps;  

        Deriv_UB_GMM      = invDeriv_eps + SE_EQ_TCC*percentile_se(2);
        Deriv_LB_GMM      = invDeriv_eps - SE_EQ_TCC*percentile_se(2);
        size(Deriv_UB_GMM)

        % Plot Data - Elasticities
        p=plot(d.x,(d.sigma-1).*Deriv_LB_CC(:,:),'-',...
         d.x,(d.sigma-1).*invDeriv_eps,'-',...
         d.x,(d.sigma-1).*Deriv_UB_CC(:,:),'-',...
         d.x,(d.sigma-1).*Deriv_LB_GMM,'.',...
         d.x,(d.sigma-1).*Deriv_UB_GMM,'.');
        xlim([-9 -4])

        standard_errors.Deriv_predicted = (d.sigma-1).*invDeriv_eps;
        standard_errors.Deriv_LB_CC     = (d.sigma-1).*Deriv_LB_CC;
        standard_errors.Deriv_UB_CC     = (d.sigma-1).*Deriv_UB_CC;
        standard_errors.Deriv_LB_GMM    = (d.sigma-1).*Deriv_LB_GMM;
        standard_errors.Deriv_UB_GMM    = (d.sigma-1).*Deriv_UB_GMM;
        
        

        
end

