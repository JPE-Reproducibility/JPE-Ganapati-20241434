function [ epsilon_elast,rho_elast,theta,x ] = make_elasticity_combo(x,k,d,o,sigma,interaction)

    switch interaction
        case 'wealth_oI'
            knots = size(k,2)/2;
            x= d.R_nE_ij;
            
            Deriv = zeros(size(x,1),2);
            [Xknots1a,k1a,Deriv1a] = evknots(knots,x(d.category==0),'cubic') ;
            [Xknots1b,k1b,Deriv1b] = evknots(knots,x(d.category==1),'cubic') ;
            d.k = [k1a k1b];
            Xknots(d.category==0,:)=Xknots1a;
            Xknots(d.category==1,:)=Xknots1b;
            Deriv(d.category==0,:)=Deriv1a;
            Deriv(d.category==1,:)=Deriv1b;
            dc0 = repmat(d.category==0,1,size(Xknots,2));
            dc1 = repmat(d.category==1,1,size(Xknots,2));
            Deriv = [Deriv.*dc0 Deriv.*dc1];

        case 'wealth_oI0'
            
            [~,~,DerivA]=evknots_make(d.k1a,x);
            Deriv = [DerivA DerivA*0];

        case 'wealth_oI1'

            [~,~,DerivA]=evknots_make(d.k1b,x);
            Deriv = [DerivA*0 DerivA];

        case 'wealth_odI'
            knots = size(k,2)/4;
            x= d.R_nE_ij;
            Deriv = zeros(size(x,1),2);
            [Xknots1a,k1a,Deriv1a] = evknots(knots,x(d.category==0),'cubic') ;
            [Xknots1b,k1b,Deriv1b] = evknots(knots,x(d.category==1),'cubic') ;
            [Xknots1c,k1c,Deriv1c] = evknots(knots,x(d.category==2),'cubic') ;
            [Xknots1d,k1d,Deriv1d] = evknots(knots,x(d.category==3),'cubic') ;
            d.k = [k1a k1b k1c k1d];
            Xknots(d.category==0,:)=Xknots1a;
            Xknots(d.category==1,:)=Xknots1b;
            Xknots(d.category==2,:)=Xknots1c;
            Xknots(d.category==3,:)=Xknots1d;

            Deriv(d.category==0,:)=Deriv1a;
            Deriv(d.category==1,:)=Deriv1b;
            Deriv(d.category==2,:)=Deriv1c;
            Deriv(d.category==3,:)=Deriv1d;
%             Deriv = [Deriv  Deriv.*repmat(d.rich_I,1,size(d.Xknots,2)) Deriv.*repmat(d.rich_J,1,size(d.Xknots,2)) Deriv.*repmat(d.rich_I.*d.rich_J,1,size(d.Xknots,2))];

            dc0 = repmat(d.category==0,1,size(Xknots,2));
            dc1 = repmat(d.category==1,1,size(Xknots,2));
            dc2 = repmat(d.category==2,1,size(Xknots,2));
            dc3 = repmat(d.category==3,1,size(Xknots,2));
            Deriv = [Deriv.*dc0 Deriv.*dc1 Deriv.*dc2 Deriv.*dc3];

        otherwise
            [~,~,Deriv]=evknots_make(k,x);
    end


    R_eps=1:size(Deriv,2);
    R_rho=size(Deriv,2)+1:2*size(Deriv,2);

    o.PRED_eps = Deriv*o.ests2(R_eps)';
    o.PRED_rho = Deriv*o.ests2(R_rho)';

    o.SE_EQ_eps = zeros(size(Deriv,1),1);
    o.SE_EQ_rho = zeros(size(Deriv,1),1);

    for i=1:size(Deriv,1)
        o.SE_EQ_eps(i) = (Deriv(i,:)*6*o.V(R_eps,R_eps)*Deriv(i,:)')^.5;
        o.SE_EQ_rho(i) = (Deriv(i,:)*6*o.V(R_rho,R_rho)*Deriv(i,:)')^.5;
    end

    percentile_se = norminv([0.025 0.975]);
    
    o.UB_eps = o.PRED_eps + o.SE_EQ_eps*percentile_se(2);
    o.LB_eps = o.PRED_eps - o.SE_EQ_eps*percentile_se(2);
    
    o.UB_rho = o.PRED_rho + o.SE_EQ_rho*percentile_se(2);
    o.LB_rho = o.PRED_rho - o.SE_EQ_rho*percentile_se(2);
    
    epsilon_elast = o.PRED_eps;
    rho_elast = o.PRED_rho;         
    theta = (1-sigma)*(1+rho_elast - epsilon_elast)./epsilon_elast;

end
            
            
                
            
            