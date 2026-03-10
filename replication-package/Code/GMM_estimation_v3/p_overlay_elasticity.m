function [ a ] = p_overlay_elasticity( d,o,type,filename,gtitle )
%P_ALPHA Summary of this function goes here
%   Detailed explanation goes here
% Create SE of predictions

%% Colors
    if filename == "../../Output/F3d5_4origin_"
        C1 = [.5 .25  0];         
        C3 = [.5 .5 1];           
        C2 = max(min(1,C1-.1),0);
        C4 = max(min(1,C3-.1),0);
    elseif filename == "../../Output/F3d1_Cross_"
        C3 = [.5 .25  0];         
        C1 = [.5 .5 1];           
        C2 = max(min(1,C1-.1),0);
        C4 = max(min(1,C3-.1),0);
    else
        C1 = [.5 .25  0];         
        C3 = [.5 .5 1];           
        C2 = max(min(1,C1-.1),0);
        C4 = max(min(1,C3-.1),0);
    end

%% Set up Matrices
    R1=1:size(d.Deriv,2);
    R2=size(d.Deriv,2)+1:2*size(d.Deriv,2);
    if nargin <= 4
        gtitle = "";
    end
    
    o.epsilon   = d.Deriv*o.ests2(R1)';
    o.rho       = d.Deriv*o.ests2(R2)';
    epsilon     = o.epsilon;
    rho         = o.rho;
    o.theta     = (d.sigma-1)*(1-(1+rho)./epsilon);

%% Standard Errors
    % percentile_se = norminv([0.025 0.975]);
    percentile_se = norminv([0.05 0.95]);
    
    o.SE_EQ_1   = zeros(size(d.Deriv,1),1);
    o.SE_EQ_2   = zeros(size(d.Deriv,1),1);
    o.SE_EQ_T   = zeros(size(d.Deriv,1),1);
    
    for i=1:size(d.Deriv,1)
        V = o.V([R1 R2],[R1 R2]);
        o.SE_EQ_1(i) = (d.Deriv(i,:)*o.V(R1,R1)*d.Deriv(i,:)')^.5;
        o.SE_EQ_2(i) = (d.Deriv(i,:)*o.V(R2,R2)*d.Deriv(i,:)')^.5;
        grad = [(1-d.sigma)*(-1/epsilon(i)^2-rho(i)/epsilon(i)^2) ; (1-d.sigma)/epsilon(i)];
        o.SE_EQ_T(i) = (grad'*(repmat(d.Deriv(i,:),2,2)*V*repmat(d.Deriv(i,:),2,2)')*grad)^.5;
    end
    
    
    
    o.UB_eps = o.epsilon + o.SE_EQ_1*percentile_se(2);
    o.LB_eps = o.epsilon - o.SE_EQ_1*percentile_se(2);
    
    o.UB_rho = o.rho + o.SE_EQ_2*percentile_se(2);
    o.LB_rho = o.rho - o.SE_EQ_2*percentile_se(2);
    
    o.UB_T = o.theta + o.SE_EQ_T*percentile_se(2);
    o.LB_T = o.theta - o.SE_EQ_T*percentile_se(2);
    
    % Extensive Margin Delta Method
    for j=1:size(d.Deriv,1)
        grad = (1/epsilon(j).^2) ;
        o.SE_EQ_T(j) = (grad'*(repmat(d.Deriv(j,:),1,1)*o.V(R1 ,R1 )*repmat(d.Deriv(j,:),1,1)')*grad)^.5;
    end
    o.extensive    = -1./o.epsilon;
    o.UB_extensive = -1./o.epsilon + o.SE_EQ_T*percentile_se(2);
    o.LB_extensive = -1./o.epsilon - o.SE_EQ_T*percentile_se(2);

    % Intensive Margin Delta Method
    for j=1:size(d.Deriv,1)
        grad = [(rho(j)/epsilon(j).^2) ; (-1/epsilon(j))];
        o.SE_EQ_T(j) = (grad'*(repmat(d.Deriv(j,:),2,2)*o.V([R1 R2],[R1 R2])*repmat(d.Deriv(j,:),2,2)')*grad)^.5;
    end
    o.UB_intensive = 1-o.rho./o.epsilon + o.SE_EQ_T*percentile_se(2);
    o.LB_intensive = 1-o.rho./o.epsilon - o.SE_EQ_T*percentile_se(2);

%% Order Estimates
    [~,b]=sort(d.R_nE_ij);
    
    d.R_nE_ij = d.R_nE_ij(b);
    o.R_nE_ij = d.R_nE_ij;
    
    try
        d.category = d.category(b);
        o.category = d.category;
    catch
    end
    
    o.theta = o.theta(b);
    o.LB_T = o.LB_T(b);
    o.UB_T = o.UB_T(b);
    
    o.epsilon = o.epsilon(b);
    o.LB_eps = o.LB_eps(b);
    o.UB_eps = o.UB_eps(b);
    
    o.rho = o.rho(b);
    o.LB_rho = o.LB_rho(b);
    o.UB_rho = o.UB_rho(b);
    
    o.extensive    = (d.sigma-1).*o.extensive(b);
    o.UB_extensive = (d.sigma-1).*o.UB_extensive(b);
    o.LB_extensive = (d.sigma-1).*o.LB_extensive(b);

    o.intensive    = (1-d.sigma)*(1-o.rho./o.epsilon);
    o.UB_intensive = (1-d.sigma)*o.UB_intensive(b);
    o.LB_intensive = (1-d.sigma)*o.LB_intensive(b);

switch type
    case 'epsilon'
        a = plot(d.R_nE_ij,o.UB_eps,'--','LineWidth',1,'Color',[.5 .5  .5]);
        hold on
        b2 = plot(d.R_nE_ij,o.epsilon,'-','LineWidth',3,'Color','black');
        plot(d.R_nE_ij,o.LB_eps,'--','LineWidth',1,'Color',[.5 .5  .5]);
        ylabel('Elasticity of Epsilon','FontSize',18)
        ylabel('\partial ln \epsilon (n)/ \partial ln n','FontSize',18)
        legend('off')
    case 'rho'
        a = plot(d.R_nE_ij,o.UB_rho,'--','LineWidth',1,'Color',[.5 .5  .5]);
        hold on
        b2 = plot(d.R_nE_ij,o.rho,'-','LineWidth',3,'Color','black');
        plot(d.R_nE_ij,o.LB_rho,'--','LineWidth',1,'Color',[.5 .5  .5]);
%         ylabel('Elasticity of Rho','FontSize',18)
%         ylabel('\partial ln \rho (n)/ \partial ln n','FontSize',18)
        ylabel('$\theta^{c}$','FontSize',24,'Interpreter','Latex')

    case 'epsilon_CC'
        a = plot(d.x,o.standard_errors1.Deriv_UB_CC,':','LineWidth',1,'Color',[.5 .5  .5]);
        hold on
        b2 = plot(d.x,o.standard_errors1.Deriv_predicted,'-','LineWidth',3,'Color','black');
        plot(d.x,o.standard_errors1.Deriv_LB_CC,':','LineWidth',1,'Color',[.5 .5  .5]);
        ylabel('Elasticity of Epsilon','FontSize',18)
        ylabel('\partial ln \epsilon (n)/ \partial ln n','FontSize',18)
        legend('off')
        b3=plot(d.R_nE_ij,o.UB_eps,'--','LineWidth',1,'Color',[.5 .5  .5]);
        plot(d.R_nE_ij,o.LB_eps,'--','LineWidth',1,'Color',[.5 .5  .5]);
        legend([ a b3 ],{'CC SE','GMM SE'},'FontSize',18,'NumColumns',1,'Location','best')
    case 'rho_CC'
        a = plot(d.x,o.standard_errors2.Deriv_UB_CC,':','LineWidth',1,'Color',[.5 .5  .5]);
        hold on
        b2 = plot(d.x,o.standard_errors2.Deriv_predicted,'-','LineWidth',3,'Color','black');
        plot(d.x,o.standard_errors2.Deriv_LB_CC,':','LineWidth',1,'Color',[.5 .5  .5]);
        ylabel('Elasticity of Rho','FontSize',18)
        ylabel('$\theta^{c}$','FontSize',24,'Interpreter','Latex')

        legend('off')
        b3=plot(d.R_nE_ij,o.UB_rho,'--','LineWidth',1,'Color',[.5 .5  .5]);
        plot(d.R_nE_ij,o.LB_rho,'--','LineWidth',1,'Color',[.5 .5  .5]);
        legend([ a b3 ],{'CC SE','GMM SE'},'FontSize',18,'NumColumns',1,'Location','best')

        legend('boxoff','Location', 'Best')
    case 'intensive'
        a = plot(d.R_nE_ij,o.UB_intensive,'--','LineWidth',1,'Color',[.5 .5  .5]);
        hold on
        b2 = plot(d.R_nE_ij,o.intensive,'-','LineWidth',3,'Color','black');
        plot(d.R_nE_ij,o.LB_intensive,'--','LineWidth',1,'Color',[.5 .5  .5]);

        ylabel('Intensive Margin','FontSize',18)
        ylabel('$d \ln \bar{x}_{ij} / d \ln \tau_{ij}$','FontSize',18,'Interpreter','Latex')

    case 'extensive'

        a = plot(d.R_nE_ij,o.UB_extensive,'--','LineWidth',1,'Color',[.5 .5  .5]);
        hold on
        b2 = plot(d.R_nE_ij,o.extensive,'-','LineWidth',3,'Color','black');
        plot(d.R_nE_ij,o.LB_extensive,'--','LineWidth',1,'Color',[.5 .5  .5]);

%         ylabel('Extensive Margin','FontSize',18)
%         ylabel('$d \ln n_{ij} / d \ln \tau_{ij}$','FontSize',18,'Interpreter','Latex')
        ylabel('$\theta^{e}(\sigma-1)$','FontSize',24,'Interpreter','Latex')
    case 'extensive_CC'
        a = plot(d.x,o.standard_errors_theta_e.Deriv_UB_CC,':','LineWidth',1,'Color',[.5 .5  .5]);
        hold on
        b2 = plot(d.R_nE_ij,o.extensive,'-','LineWidth',3,'Color','black');

%         b2 = plot(d.x,o.standard_errors_theta_e.Deriv_predicted,'-','LineWidth',3,'Color','black');
        plot(d.x,o.standard_errors_theta_e.Deriv_LB_CC,':','LineWidth',1,'Color',[.5 .5  .5]);
        ylabel('$\theta^{e}(\sigma-1)$','FontSize',24,'Interpreter','Latex')
        legend('off')
        b3=plot(d.R_nE_ij,o.UB_extensive,'--','LineWidth',1,'Color',[.5 .5  .5]);
        plot(d.R_nE_ij,o.LB_extensive,'--','LineWidth',1,'Color',[.5 .5  .5]);
        legend([ a b3 ],{'CC SE','GMM SE'},'FontSize',18,'NumColumns',1,'Location','best')

    case 'epsilon_split'
        a = plot(d.R_nE_ij(d.category==1),o.UB_eps(d.category==1),'--','LineWidth',1,'Color',C3);
        hold on
        b2 = plot(d.R_nE_ij(d.category==1),o.epsilon(d.category==1),'-','LineWidth',3,'Color',C3);
        plot(d.R_nE_ij(d.category==1),o.LB_eps(d.category==1),'--','LineWidth',1,'Color',C3);
        plot(d.R_nE_ij(d.category==0),o.UB_eps(d.category==0),'--','LineWidth',1,'Color',C1);
        b5 = plot(d.R_nE_ij(d.category==0),o.epsilon(d.category==0),'-','LineWidth',3,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_eps(d.category==0),'--','LineWidth',1,'Color',C1);
        ylabel('Elasticity of Epsilon','FontSize',18)
        ylabel('\partial ln \epsilon (n)/ \partial ln n','FontSize',18)

        legend([ b5 b2],{[d.type0],[d.type1]},'FontSize',18,'NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'rho_split'
        a = plot(d.R_nE_ij(d.category==1),o.UB_rho(d.category==1),'--','LineWidth',1,'Color',C3);
        hold on
        b2 = plot(d.R_nE_ij(d.category==1),o.rho(d.category==1),'-','LineWidth',3,'Color',C3);
        plot(d.R_nE_ij(d.category==1),o.LB_rho(d.category==1),'--','LineWidth',1,'Color',C3);
        plot(d.R_nE_ij(d.category==0),o.UB_rho(d.category==0),'--','LineWidth',1,'Color',C1);
        b5 = plot(d.R_nE_ij(d.category==0),o.rho(d.category==0),'-','LineWidth',3,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_rho(d.category==0),'--','LineWidth',1,'Color',C1);
%         ylabel('Elasticity of Rho','FontSize',18)
%         ylabel('\partial ln \rho (n)/ \partial ln n','FontSize',18)
        ylabel('$\theta^{c}$','FontSize',24,'Interpreter','Latex')

        legend([ b5 b2],{[d.type0],[d.type1]},'FontSize',18,'NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'intensive_split'
        a = plot(d.R_nE_ij(d.category==1),o.UB_intensive(d.category==1),'--','LineWidth',1,'Color',C3);
        hold on
        b2 = plot(d.R_nE_ij(d.category==1),o.intensive(d.category==1),'-','LineWidth',3,'Color',C3);
        plot(d.R_nE_ij(d.category==1),o.LB_intensive(d.category==1),'--','LineWidth',1,'Color',C3);
        plot(d.R_nE_ij(d.category==0),o.UB_intensive(d.category==0),'--','LineWidth',1,'Color',C1);
        b5 = plot(d.R_nE_ij(d.category==0),o.intensive(d.category==0),'-','LineWidth',3,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_intensive(d.category==0),'--','LineWidth',1,'Color',C1);
        ylabel('Intensive Margin','FontSize',18)
        ylabel('$d \ln \bar{x}_{ij} / d \ln \tau_{ij}$','FontSize',18,'Interpreter','Latex')

        legend([ b5 b2],{[d.type0],[d.type1]},'FontSize',18,'NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'extensive_split'

        a = plot(d.R_nE_ij(d.category==1),o.UB_extensive(d.category==1),'--','LineWidth',1,'Color',C3);
        hold on
        b2 = plot(d.R_nE_ij(d.category==1),o.extensive(d.category==1),'-','LineWidth',3,'Color',C3);
        plot(d.R_nE_ij(d.category==1),o.LB_extensive(d.category==1),'--','LineWidth',1,'Color',C3);
        plot(d.R_nE_ij(d.category==0),o.UB_extensive(d.category==0),'--','LineWidth',1,'Color',C1);
        b5 = plot(d.R_nE_ij(d.category==0),o.extensive(d.category==0),'-','LineWidth',3,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_extensive(d.category==0),'--','LineWidth',1,'Color',C1);
%         ylabel('Extensive Margin','FontSize',18)
%         ylabel('$d \ln n_{ij} / d \ln \tau_{ij}$','FontSize',18,'Interpreter','Latex')
        ylabel('$\theta^{e}(\sigma-1)$','FontSize',24,'Interpreter','Latex')

        legend([ b5 b2],{[d.type0],[d.type1]},'FontSize',18,'NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'theta'
        a = plot(d.R_nE_ij,o.theta,'-','LineWidth',3,'Color','black');
        hold on
        b2 = plot(d.R_nE_ij,o.UB_T,'--','LineWidth',1,'Color',[.5 .5  .5]);
        plot(d.R_nE_ij,o.LB_T,'--','LineWidth',1,'Color',[.5 .5  .5]);
        hold off
%         ylabel('-\partial ln X_{ij} / \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')

        ylim([min([1.5 min(o.theta)]) max([10 max(o.theta)])]);

    case 'theta_split'
        a = plot(d.R_nE_ij(d.category==0),o.theta(d.category==0),'-','LineWidth',3,'Color',C1);
        hold on
        plot(d.R_nE_ij(d.category==0),o.UB_T(d.category==0),'--','LineWidth',1,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_T(d.category==0),'--','LineWidth',1,'Color',C1);
        plot(d.R_nE_ij(d.category==1),o.UB_T(d.category==1),'--','LineWidth',1,'Color',C3);
        plot(d.R_nE_ij(d.category==1),o.LB_T(d.category==1),'--','LineWidth',1,'Color',C3);
        b5 = plot(d.R_nE_ij(d.category==1),o.theta(d.category==1),'-','LineWidth',3,'Color',C3);
%         ylabel('-\partial ln X_{ij} / \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')
       
        legend([a b5 ],{[d.type0],[d.type1]},'FontSize',18,'NumColumns',1,'Location','best')
        ylim([min([1.5 min(o.theta)]) max([10 max(o.theta)])]);
        legend('boxoff','Location', 'Best')

    case 'rho_quad'
        a = plot(d.R_nE_ij(d.category==0),o.UB_rho(d.category==0),'--','LineWidth',.5,'Color',C1);
        hold on
        b2 = plot(d.R_nE_ij(d.category==0),o.rho(d.category==0),'-','LineWidth',2,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_rho(d.category==0),'--','LineWidth',.5,'Color',C1);
        plot(d.R_nE_ij(d.category==1),o.UB_rho(d.category==1),'--','LineWidth',.5,'Color',C2);
        b5 = plot(d.R_nE_ij(d.category==1),o.rho(d.category==1),'-','LineWidth',4.5,'Color',C2);
        plot(d.R_nE_ij(d.category==1),o.LB_rho(d.category==1),'--','LineWidth',.5,'Color',C2);
        plot(d.R_nE_ij(d.category==2),o.UB_rho(d.category==2),'--','LineWidth',.5,'Color',C3);
        b6 = plot(d.R_nE_ij(d.category==2),o.rho(d.category==2),'-','LineWidth',2,'Color',C3);
        plot(d.R_nE_ij(d.category==2),o.LB_rho(d.category==2),'--','LineWidth',.5,'Color',C3);
        plot(d.R_nE_ij(d.category==3),o.UB_rho(d.category==3),'--','LineWidth',.5,'Color',C4);
        b7 = plot(d.R_nE_ij(d.category==3),o.rho(d.category==3),'-','LineWidth',4.5,'Color',C4);
        plot(d.R_nE_ij(d.category==3),o.LB_rho(d.category==3),'--','LineWidth',.5,'Color',C4);
%         ylabel('\partial ln \rho (n)/ \partial ln n','FontSize',18)
        ylabel('$\theta^{c}$','FontSize',24,'Interpreter','Latex')
        legend([b6 b7 b2 b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'epsilon_quad'
        a = plot(d.R_nE_ij(d.category==0),o.UB_eps(d.category==0),'--','LineWidth',.5,'Color',C1);
        hold on
        b2 = plot(d.R_nE_ij(d.category==0),o.epsilon(d.category==0),'-','LineWidth',2,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_eps(d.category==0),'--','LineWidth',.5,'Color',C1);
        plot(d.R_nE_ij(d.category==1),o.UB_eps(d.category==1),'--','LineWidth',.5,'Color',C2);
        b5 = plot(d.R_nE_ij(d.category==1),o.epsilon(d.category==1),'-','LineWidth',4.5,'Color',C2);
        plot(d.R_nE_ij(d.category==1),o.LB_eps(d.category==1),'--','LineWidth',.5,'Color',C2);
        plot(d.R_nE_ij(d.category==2),o.UB_eps(d.category==2),'--','LineWidth',.5,'Color',C3);
        b6 = plot(d.R_nE_ij(d.category==2),o.epsilon(d.category==2),'-','LineWidth',2,'Color',C3);
        plot(d.R_nE_ij(d.category==2),o.LB_eps(d.category==2),'--','LineWidth',.5,'Color',C3);
        plot(d.R_nE_ij(d.category==3),o.UB_eps(d.category==3),'--','LineWidth',.5,'Color',C4);
        b7 = plot(d.R_nE_ij(d.category==3),o.epsilon(d.category==3),'-','LineWidth',4.5,'Color',C4);
        plot(d.R_nE_ij(d.category==3),o.LB_eps(d.category==3),'--','LineWidth',.5,'Color',C4);
        ylabel('\partial ln \epsilon (n) / \partial ln n','FontSize',18)

        legend([b6 b7 b2 b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'theta_quad'
        b2 = plot(d.R_nE_ij(d.category==0),o.theta(d.category==0),'-','LineWidth',2,'Color',C1);
        hold on
        a = plot(d.R_nE_ij(d.category==0),o.UB_T(d.category==0),'--','LineWidth',.5,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_T(d.category==0),'--','LineWidth',.5,'Color',C1);
        plot(d.R_nE_ij(d.category==1),o.UB_T(d.category==1),'--','LineWidth',.5,'Color',C2);
        b5 = plot(d.R_nE_ij(d.category==1),o.theta(d.category==1),'-','LineWidth',4.5,'Color',C2);
        plot(d.R_nE_ij(d.category==1),o.LB_T(d.category==1),'--','LineWidth',.5,'Color',C2);
        plot(d.R_nE_ij(d.category==2),o.UB_T(d.category==2),'--','LineWidth',.5,'Color',C3);
        b6 = plot(d.R_nE_ij(d.category==2),o.theta(d.category==2),'-','LineWidth',2,'Color',C3);
        plot(d.R_nE_ij(d.category==2),o.LB_T(d.category==2),'--','LineWidth',.5,'Color',C3);
        plot(d.R_nE_ij(d.category==3),o.UB_T(d.category==3),'--','LineWidth',.5,'Color',C4);
        b7 = plot(d.R_nE_ij(d.category==3),o.theta(d.category==3),'-','LineWidth',4.5,'Color',C4);
        plot(d.R_nE_ij(d.category==3),o.LB_T(d.category==3),'--','LineWidth',.5,'Color',C4);
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')

        if filename == "output/F3d5_4origin_quad_"
            legend([b2 b5 b6 b7 ],{[d.type0],[d.type1],[d.type2],[d.type3]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
            legend('boxoff','Location', 'southwest')
        else
            legend([b6 b7 b2 b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
            legend('boxoff','Location', 'Best')
        end

        ylim([min([0 min(o.theta)]) max([10 max(o.theta)])]);
    case 'theta_quadA'
        plot(d.R_nE_ij(d.category==2),o.UB_T(d.category==2),'--','LineWidth',.5,'Color',C3);
        hold on
        a = plot(d.R_nE_ij(d.category==2),o.theta(d.category==2),'-','LineWidth',2,'Color',C3);
        plot(d.R_nE_ij(d.category==2),o.LB_T(d.category==2),'--','LineWidth',.5,'Color',C3);
        plot(d.R_nE_ij(d.category==3),o.UB_T(d.category==3),'--','LineWidth',.5,'Color',C4);
        b7 = plot(d.R_nE_ij(d.category==3),o.theta(d.category==3),'-','LineWidth',4.5,'Color',C4);
        plot(d.R_nE_ij(d.category==3),o.LB_T(d.category==3),'--','LineWidth',.5,'Color',C4);
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')
        legend([a b7   ],{[d.type2],[d.type3]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')
        ylim([min([0 min(o.theta)]) max([10 max(o.theta)])]);
    case 'theta_quadB'
        b2 = plot(d.R_nE_ij(d.category==0),o.theta(d.category==0),'-','LineWidth',2,'Color',C1);
        hold on
        a = plot(d.R_nE_ij(d.category==0),o.UB_T(d.category==0),'--','LineWidth',.5,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_T(d.category==0),'--','LineWidth',.5,'Color',C1);
        plot(d.R_nE_ij(d.category==1),o.UB_T(d.category==1),'--','LineWidth',.5,'Color',C2);
        b5 = plot(d.R_nE_ij(d.category==1),o.theta(d.category==1),'-','LineWidth',4.5,'Color',C2);
        plot(d.R_nE_ij(d.category==1),o.LB_T(d.category==1),'--','LineWidth',.5,'Color',C2);
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')
        legend([  b2 b5 ],{[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')
        ylim([min([0 min(o.theta)]) max([10 max(o.theta)])]);
    case 'theta_quad_nose'
        a = plot(d.R_nE_ij(d.category==0),o.theta(d.category==0),'-','LineWidth',2,'Color',C1);
        hold on
        b5 = plot(d.R_nE_ij(d.category==1),o.theta(d.category==1),'-','LineWidth',4.5,'Color',C2);
        b6 = plot(d.R_nE_ij(d.category==2),o.theta(d.category==2),'-','LineWidth',2,'Color',C3);
        b7 = plot(d.R_nE_ij(d.category==3),o.theta(d.category==3),'-','LineWidth',4.5,'Color',C4);
%         ylabel('-\partial ln X_{ij} / \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')

        legend([b6 b7 a b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

        ylim([min([0 min(o.theta)]) max([10 max(o.theta)])]);

    case 'intensive_quad'
        a = plot(d.R_nE_ij(d.category==0),o.UB_intensive(d.category==0),'--','LineWidth',1,'Color',C1);
        hold on
        b2 = plot(d.R_nE_ij(d.category==0),o.intensive(d.category==0),'--','LineWidth',3.5,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_intensive(d.category==0),'--','LineWidth',1,'Color',C1);
        plot(d.R_nE_ij(d.category==1),o.UB_intensive(d.category==1),'--','LineWidth',1,'Color',C2);
        b5 = plot(d.R_nE_ij(d.category==1),o.intensive(d.category==1),':','LineWidth',3,'Color',C2);
        plot(d.R_nE_ij(d.category==1),o.LB_intensive(d.category==1),'--','LineWidth',1,'Color',C2);
        plot(d.R_nE_ij(d.category==2),o.UB_intensive(d.category==2),'--','LineWidth',1,'Color',C3);
        b6 = plot(d.R_nE_ij(d.category==2),o.intensive(d.category==2),'-','LineWidth',2.5,'Color',C3);
        plot(d.R_nE_ij(d.category==2),o.LB_intensive(d.category==2),'--','LineWidth',1,'Color',C3);
        plot(d.R_nE_ij(d.category==3),o.UB_intensive(d.category==3),'--','LineWidth',1,'Color',C4);
        b7 = plot(d.R_nE_ij(d.category==3),o.intensive(d.category==3),'-.','LineWidth',2,'Color',C4);
        plot(d.R_nE_ij(d.category==3),o.LB_intensive(d.category==3),'--','LineWidth',1,'Color',C4);
        ylabel('Intensive Margin','FontSize',18)
        ylabel('$d \ln \bar{x}_{ij} / d \ln \tau_{ij}$','FontSize',18,'Interpreter','Latex')

        legend([b6 b7 b2 b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'extensive_quad'
        a = plot(d.R_nE_ij(d.category==0),o.UB_extensive(d.category==0),'--','LineWidth',.5,'Color',C1);
        hold on
        b2 = plot(d.R_nE_ij(d.category==0),o.extensive(d.category==0),'-','LineWidth',2,'Color',C1);
        plot(d.R_nE_ij(d.category==0),o.LB_extensive(d.category==0),'--','LineWidth',.5,'Color',C1);
        plot(d.R_nE_ij(d.category==1),o.UB_extensive(d.category==1),'--','LineWidth',.5,'Color',C2);
        b5 = plot(d.R_nE_ij(d.category==1),o.extensive(d.category==1),'-','LineWidth',5,'Color',C2);
        plot(d.R_nE_ij(d.category==1),o.LB_extensive(d.category==1),'--','LineWidth',.5,'Color',C2);
        plot(d.R_nE_ij(d.category==2),o.UB_extensive(d.category==2),'--','LineWidth',.5,'Color',C3);
        b6 = plot(d.R_nE_ij(d.category==2),o.extensive(d.category==2),'-','LineWidth',2,'Color',C3);
        plot(d.R_nE_ij(d.category==2),o.LB_extensive(d.category==2),'--','LineWidth',.5,'Color',C3);
        plot(d.R_nE_ij(d.category==3),o.UB_extensive(d.category==3),'--','LineWidth',.5,'Color',C4);
        b7 = plot(d.R_nE_ij(d.category==3),o.extensive(d.category==3),'-','LineWidth',5,'Color',C4);
        plot(d.R_nE_ij(d.category==3),o.LB_extensive(d.category==3),'--','LineWidth',.5,'Color',C4);
%         ylabel('Extensive Margin','FontSize',18)
%         ylabel('\partial ln n_{ij)/ \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta^{e}(\sigma-1)$','FontSize',24,'Interpreter','Latex')

        legend([b6 b7 b2 b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')


    case 'theta_quadl'
        a = plot(d.R_nE_ij(d.category==0),o.theta(d.category==0),'--','LineWidth',3,'Color',C1);
        hold on
        b5 = plot(d.R_nE_ij(d.category==1),o.theta(d.category==1),'-','LineWidth',3,'Color',C2);
        b6 = plot(d.R_nE_ij(d.category==2),o.theta(d.category==2),'-','LineWidth',3,'Color',C3);
        b7 = plot(d.R_nE_ij(d.category==3),o.theta(d.category==3),'-.','LineWidth',3,'Color',C4);
%         ylabel('-\partial ln X_{ij} / \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')
        set(gca, 'YScale', 'log')
        legend([b6 b7 a b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')

    case 'theta_quadline'
        a = plot(d.R_nE_ij(d.category==0),o.theta(d.category==0),'--','Color',C1,'LineWidth',2);
        hold on
        b5 = plot(d.R_nE_ij(d.category==1),o.theta(d.category==1),'-','Color',C2,'LineWidth',2);
        b6 = plot(d.R_nE_ij(d.category==2),o.theta(d.category==2),'-','Color',C3,'LineWidth',2);
        b7 = plot(d.R_nE_ij(d.category==3),o.theta(d.category==3),'-.','Color',C4,'LineWidth',2);
%         ylabel('-\partial ln X_{ij} / \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')
        set(gca, 'YScale', 'log')
        legend([b6 b7 a b5 ],{[d.type2],[d.type3],[d.type0],[d.type1]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'Best')


end

%% Clean Up
    hold off
    xlabel('Log Exporter Firm Share','FontSize',18)

    ax = gca;
    ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];
    ax.XTickLabel = {'0.01%' '0.1%','1%','10%','100%'};
    ax.FontSize = 18;
    xlim([-9.9 -2]);
    
    if filename == "output/F3d5_4origin_"
        legend([b2 b5 b6 b7 ],{[d.type0],[d.type1],[d.type2],[d.type3]},'FontSize',18,'Location','northeast','NumColumns',1,'Location','best')
        legend('boxoff','Location', 'southwest')
    end
    
    title(gtitle);
    saveas(a,[filename,type],'epsc')
    
    a = o.theta;

end

