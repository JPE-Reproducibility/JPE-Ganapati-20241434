function [standard_errors] = CC_se_splines(R1,ests2,d,uhat,G_bar,W,Lambda,V)
%% Standard Errors full functional form
% Code From CC 2017 QE

%         [Psixx,~,Deriv] =evknots_make(d.k,d.x);
    
        % Cubic B-splines, plot over 5%-95% empirical quantiles
        r = d.r;
        xtemp = (0.05:0.001:0.95)';
%         xtemp = (0.01:0.001:0.99)';

        [Psixx, Deriv] = bspline(xtemp, 0, r, d.k);
        % Rescale because of normalization to [0, 1] as before, trimming
        % very extreme quantiles
%         scale_upper = quantile(d.R_nE_ij, 0.999);
%         scale_lower = quantile(d.R_nE_ij, 0.001);

        scale_upper = d.scale_upper;
        scale_lower = d.scale_lower;
%         scale_upper = 1; %%%%%%%%
%         scale_lower = 0; %%%%%%%%

        d.x = xtemp * (scale_upper - scale_lower) + scale_lower;
        % Rescale derivative 
        Deriv = Deriv / (scale_upper - scale_lower);
        
        Nboot = 1000;
        SGl   = pinv(G_bar(:,R1)'*W*G_bar(:,R1))*G_bar(:,R1)'*W;
        Q     = SGl*Lambda*SGl';
        nxx   = length(Psixx);
        Vxx   = zeros(nxx,1);
        Sxx   = zeros(nxx,1);
        B = d.ZA_K;
        [n,K] = size(B);
        alpha = [ 0.9500 ];
                
        for i = 1:nxx
            Vxx(i) = Psixx(i,:)*Q*Psixx(i,:)';
            Sxx(i) = sqrt(Vxx(i));
        end
     
        SGl   = SGl(:,1:size(SGl,2)/2);
        
        %   bootsrap subroutine
        zb   = zeros(Nboot,1);
        for i = 1:Nboot
            % Mammen two-point
            w     = 0.5 - 0.5*sqrt(5)*(2*(rand(n,1) < (sqrt(5)+1)/(2*sqrt(5)))-1);
            R     = SGl*(B'*(uhat.*w))/sqrt(n);
            Rb    = (Psixx*R)./Sxx;
            zb(i) = max(abs(Rb));
        end
        zalpha = quantile(zb,alpha);

        c       = ests2(R1);
        Prediction  = Psixx*c';
        v       = Sxx/sqrt(n)*zalpha;
        hxxal   = Prediction*ones(1,length(alpha));
        LB_CC   = hxxal-v;
        UB_CC   = hxxal+v;   
        
        %% For Derivative
        for i = 1:nxx
            Vxx(i) = Deriv(i,:)*Q*Deriv(i,:)';
            Sxx(i) = sqrt(6*Vxx(i));
        end
        
        zb   = zeros(Nboot,1);
        for i = 1:Nboot
            % Mammen two-point
            w     = 0.5 - 0.5*sqrt(5)*(2*(rand(n,1) < (sqrt(5)+1)/(2*sqrt(5)))-1);
            R     = SGl*(B'*(uhat.*w))/sqrt(n);
            Rb    = (Deriv*R)./Sxx;
            zb(i) = max(abs(Rb));
        end
        zalpha = quantile(zb,alpha);
        v     = Sxx/sqrt(n)*zalpha;

        deriv_hxxal = Deriv*ests2(R1)'*ones(1,length(alpha));
        Deriv_LB_CC    = deriv_hxxal-v;
        Deriv_UB_CC    = deriv_hxxal+v;                
        
        %% Create SE of predictions - GMM
        for i=1:size(Psixx,1)
            SE_EQ_1_GMM(i) = (Psixx(i,:)*1*V(R1,R1)*Psixx(i,:)')^.5;
        end
        
        percentile_se = norminv([0.025 0.975]);
        UB_GMM      = Prediction + SE_EQ_1_GMM'*percentile_se(2);
        LB_GMM      = Prediction - SE_EQ_1_GMM'*percentile_se(2);
        
%         % Plot Data - Levels
%         p=plot(d.x,UB_CC(:,:),'-',...
%                 d.x,Prediction,'.',...
%                 d.x,LB_CC(:,:),'-',...
%                 d.x,UB_GMM,'.',...
%                 d.x,Prediction,'.',...
%                 d.x,LB_GMM,'.');

        Deriv_predicted = Deriv*ests2(R1)';

        SE_EQ_Deriv_GMM = zeros(size(Deriv,1),1);

        for i=1:size(Deriv,1)
            SE_EQ_Deriv_GMM(i) = (Deriv(i,:)*6*V(R1,R1)*Deriv(i,:)')^.5;
        end
        Deriv_UB_GMM      = Deriv_predicted + SE_EQ_Deriv_GMM*percentile_se(2);
        Deriv_LB_GMM      = Deriv_predicted - SE_EQ_Deriv_GMM*percentile_se(2);
        
%         % Plot Data - Elasticities
%         p=plot(d.x,Deriv_LB_CC(:,:),'-',...
%          d.x,Deriv_predicted,'-',...
%          d.x,Deriv_UB_CC(:,:),'-',...
%          d.x,Deriv_LB_GMM,'.',...
%          d.x,Deriv_predicted,'.',...
%          d.x,Deriv_UB_GMM,'.');
     
        standard_errors.Prediction      = Prediction;
        standard_errors.LB_CC           = LB_CC;
        standard_errors.UB_CC           = UB_CC;
        standard_errors.LB_GMM          = LB_GMM;
        standard_errors.UB_GMM          = UB_GMM;
        
        standard_errors.Deriv_predicted = Deriv_predicted;
        standard_errors.Deriv_LB_CC     = Deriv_LB_CC;
        standard_errors.Deriv_UB_CC     = Deriv_UB_CC;
        standard_errors.Deriv_LB_GMM    = Deriv_LB_GMM;
        standard_errors.Deriv_UB_GMM    = Deriv_UB_GMM;

        
end

