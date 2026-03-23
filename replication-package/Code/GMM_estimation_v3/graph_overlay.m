function graph_overlay( o,d,type,filename,gtitle,lit_overlay)
%P_epsilon Summary of this function goes here
%   Detailed explanation goes here
    % Create SE of predictions
    num = size(o,2);
    for i= 1:num
        d1R1 = 1:size(d(i).X_K,2);
        d1R2=size(d(i).X_K,2)+1:2*size(d(i).X_K,2);

        o(i).epsilon = d(i).Deriv*o(i).ests2(d1R1)';
        o(i).rho = d(i).Deriv*o(i).ests2(d1R2)';

        o(i).SE_EQ_eps = zeros(size(d(i).Deriv,1),1);
        o(i).SE_EQ_rho = zeros(size(d(i).Deriv,1),1);
    
        for j=1:size(d(i).Deriv,1)
            o(i).SE_EQ_eps(j) = (d(i).Deriv(j,:)*6*o(i).V(d1R1,d1R1)*d(i).Deriv(j,:)')^.5;
            o(i).SE_EQ_rho(j) = (d(i).Deriv(j,:)*6*o(i).V(d1R2,d1R2)*d(i).Deriv(j,:)')^.5;
        end

        percentile_se = norminv([0.05 0.95]);
        
        o(i).UB_eps = o(i).epsilon + o(i).SE_EQ_eps*percentile_se(2);
        o(i).LB_eps = o(i).epsilon - o(i).SE_EQ_eps*percentile_se(2);
        
        o(i).UB_rho = o(i).rho + o(i).SE_EQ_rho*percentile_se(2);
        o(i).LB_rho = o(i).rho - o(i).SE_EQ_rho*percentile_se(2);

        rho1 = o(i).rho;         epsilon1 = o(i).epsilon;
        o(i).theta = (1-d(i).sigma)*(1+rho1 - epsilon1)./epsilon1;


        % Theta Delta Method
        for j=1:size(d(i).Deriv,1)
            o(i).SE_EQ_1(j) = (d(i).Deriv(j,:)*o(i).V(d1R1,d1R1)*d(i).Deriv(j,:)')^.5;
            o(i).SE_EQ_2(j) = (d(i).Deriv(j,:)*o(i).V(d1R2,d1R2)*d(i).Deriv(j,:)')^.5;
            grad = [(1-d(i).sigma)*(-1/epsilon1(j).^2-rho1(j)/epsilon1(j).^2) ; (1-d(i).sigma)/epsilon1(j)];
            o(i).SE_EQ_T(j) = (grad'*(repmat(d(i).Deriv(j,:),1,2)*o(i).V([d1R1 d1R2],[d1R1 d1R2])*repmat(d(i).Deriv(j,:),1,2)')*grad)^.5;
        end

        o(i).UB_T = o(i).theta + o(i).SE_EQ_T'*percentile_se(2);
        o(i).LB_T = o(i).theta - o(i).SE_EQ_T'*percentile_se(2);


        % Extensive Margin Delta Method
        for j=1:size(d(i).Deriv,1)
            grad = [(1/epsilon1(j).^2) ];
            o(i).SE_EQ_T(j) = (grad'*(repmat(d(i).Deriv(j,:),1,1)*o(i).V(d1R1 ,d1R1 )*repmat(d(i).Deriv(j,:),1,1)')*grad)^.5;
        end

        o(i).UB_extensive = -1./o(i).epsilon + o(i).SE_EQ_T'*percentile_se(2);
        o(i).LB_extensive = -1./o(i).epsilon - o(i).SE_EQ_T'*percentile_se(2);


        % Intensive Margin Delta Method
        for j=1:size(d(i).Deriv,1)
            grad = [(rho1(j)/epsilon1(j).^2) ; (-1/epsilon1(j))];
            o(i).SE_EQ_T(j) = (grad'*(repmat(d(i).Deriv(j,:),2,2)*o(i).V([d1R1 d1R2],[d1R1 d1R2])*repmat(d(i).Deriv(j,:),2,2)')*grad)^.5;
        end
        o(i).UB_intensive = 1-o(i).rho./o(i).epsilon + o(i).SE_EQ_T'*percentile_se(2);
        o(i).LB_intensive = 1-o(i).rho./o(i).epsilon - o(i).SE_EQ_T'*percentile_se(2);

    end


    z(1).color = [0 0 0];
    z(2).color = [.5 .75 .5];
    z(3).color = [.75 .5 .15];
    z(4).color = [.75 .15 .25];

    for i= 1:num
        [~,ts] = sort(d(i).R_nE_ij);
        d(i).R_nE_ij= d(i).R_nE_ij(ts);
        o(i).epsilon= o(i).epsilon(ts);
        o(i).UB_eps= o(i).UB_eps(ts);
        o(i).LB_eps= o(i).LB_eps(ts);
        o(i).rho= o(i).rho(ts);
        o(i).UB_rho= o(i).UB_rho(ts);
        o(i).LB_rho= o(i).LB_rho(ts);
        o(i).theta= o(i).theta(ts);
        o(i).UB_T= o(i).UB_T(ts);
        o(i).LB_T= o(i).LB_T(ts);
        o(i).UB_intensive= o(i).UB_intensive(ts);
        o(i).LB_intensive= o(i).LB_intensive(ts);
        o(i).UB_extensive= o(i).UB_extensive(ts);
        o(i).LB_extensive= o(i).LB_extensive(ts);

        try
        d(i).category = d(i).category(ts);
        catch
        end

    end

    
    if nargin < 6
        lit_overlay = 0;
    end
%     load('literature.mat')
%     lit_overlay = 1;



    switch type
        case 'epsilon'
            clf
            hold on
            for i= 1:num
                a(i) = plot(d(i).R_nE_ij,o(i).epsilon,'-','LineWidth',3,'Color',z(i).color);
%                  if lit_overlay == 0
                plot(d(i).R_nE_ij,o(i).UB_eps,'--','LineWidth',1,'Color',z(i).color);
                plot(d(i).R_nE_ij,o(i).LB_eps,'--','LineWidth',1,'Color',z(i).color);
%                  end
            end
            ylabel('\partial ln \epsilon (n ) / \partial ln (n)','FontSize',18)
        
            if lit_overlay == 1
                load('literature.mat')
%                 a(num+1)=plot(ln_ij,elast_e_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]	);
                a(num+1)=plot(ln_ij,elast_e_tpareto(range1),'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]	);
                a(num+2)=plot(ln_ij,elast_e_logn(range1),':','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
                a(num+3)=plot(log(unif_EKK(range2)),smoothdata(elast_e_EKK(range2),'gaussian',20),'-.','LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
            end

        case 'rho'
            clf
            hold on
            for i= 1:num
                a(i) = plot(d(i).R_nE_ij,o(i).rho,'-','LineWidth',3,'Color',z(i).color);
%                 if lit_overlay == 0
                plot(d(i).R_nE_ij,o(i).UB_rho,'--','LineWidth',1,'Color',z(i).color);
                plot(d(i).R_nE_ij,o(i).LB_rho,'--','LineWidth',1,'Color',z(i).color);
%                 end
            end

            if lit_overlay == 1
                load('literature.mat')
%                 a(num+1)=plot(ln_ij,elast_r_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]	);
                a(num+1)=plot(ln_ij,elast_r_tpareto(range1),'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]	);
                a(num+2)=plot(ln_ij,elast_r_logn(range1),':','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
                a(num+3)=plot(log(unif_EKK(range2)),smoothdata(elast_r_EKK(range2),'gaussian',20),'-.','LineWidth',2,'Color',[0.4940 0.1840 0.5560]	);
            end

%             ylabel('\partial ln \rho (n ) / \partial ln (n)','FontSize',18)
        ylabel('$\theta^{c}$','FontSize',24,'Interpreter','Latex')

        case 'theta'
            clf
            hold on
            for i= 1:num
                a(i) = plot(d(i).R_nE_ij,o(i).theta,'-','LineWidth',3,'Color',z(i).color);
%                 if lit_overlay == 0
                plot(d(i).R_nE_ij,o(i).UB_T,'--','LineWidth',1,'Color',z(i).color);
                plot(d(i).R_nE_ij,o(i).LB_T,'--','LineWidth',1,'Color',z(i).color);
%                 end
            end

            if lit_overlay == 1
                load('literature.mat')
%                 a(num+1)=plot(ln_ij,theta_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]	);
                a(num+1)=plot(ln_ij,theta_tpareto(range1),'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]	);
                a(num+2)=plot(ln_ij,theta_logn(range1),':','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
                a(num+3)=plot(log(unif_EKK(range2)),smoothdata(theta_EKK(range2),'gaussian',20),'-.','LineWidth',2,'Color',[0.4940 0.1840 0.5560]		);
            end            

        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')

        case 'theta_nose'
            clf
            hold on
            for i= 1:num
                a(i) = plot(d(i).R_nE_ij,o(i).theta,'-','LineWidth',3,'Color',z(i).color);
            end

            if lit_overlay == 1
                load('literature.mat')
%                 a(num+1)=plot(ln_ij,theta_pareto(range1),'-','LineWidth',2.5,'Color','Red');
                a(num+1)=plot(ln_ij,theta_tpareto(range1),'--','LineWidth',2.5,'Color','Green');
                a(num+2)=plot(ln_ij,theta_logn(range1),'-','LineWidth',2.5,'Color','Blue');
%                 a(num+3)=plot(ln_ij,theta_logn(range1),'-','LineWidth',2.5,'Color','Blue');
                %p4=plot(log(unif_EKK(range2)),smoothdata(theta_EKK(range2),'gaussian',10),'-','LineWidth',2.5,'Color','Orange');
            end            

%             ylabel('-\partial ln X_{ij} / \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')


    case 'theta_split'

        clf
        hold on
        for i= 1:num
            z(1).color = [.5 .25  0];         C2 = [.5 0 0];
            z(2).color = [.5 .5 1];           C4 = [0 .5 .5];
        
            C1 = min(1,z(1).color+(i-1)/1.5);
            C2 = min(1,z(2).color+(i-1)/4);


            a(i) = plot(d(i).R_nE_ij(d(i).category==0),o(i).theta(d(i).category==0),'-','LineWidth',3,'Color',C1);
            b5(i) = plot(d(i).R_nE_ij(d(i).category==1),o(i).theta(d(i).category==1),':','LineWidth',3,'Color',C2);
            plot(d(i).R_nE_ij(d(i).category==0),o(i).UB_T(d(i).category==0),'--','LineWidth',.5,'Color',C1);
            plot(d(i).R_nE_ij(d(i).category==0),o(i).LB_T(d(i).category==0),'--','LineWidth',.5,'Color',C1);
            plot(d(i).R_nE_ij(d(i).category==1),o(i).UB_T(d(i).category==1),'--','LineWidth',.5,'Color',C2);
            plot(d(i).R_nE_ij(d(i).category==1),o(i).LB_T(d(i).category==1),'--','LineWidth',.5,'Color',C2);
        end

%         ylabel('\partial ln X_{ij} / \partial ln \tau_{ij}','FontSize',18)
        ylabel('$\theta(\sigma-1)$','FontSize',24,'Interpreter','Latex')

        hold on
        legend([a(1) b5(1) a(2) b5(2)],{[d(1).type0],[d(1).type1],[d(2).type0 ' - Sectoral'],[d(2).type1 ' - Sectoral']},'FontSize',18,'NumColumns',1,'Location','best')
        ylim([min([1.5 ]) max([10 ])]);
        legend('boxoff','Location', 'Best')

    case 'intensive'
        clf
        hold on
        for i= 1:num
            a(i) = plot(d(i).R_nE_ij,(1-d(i).sigma)*(1-o(i).rho./o(i).epsilon)  ,'-','LineWidth',3,'Color',z(i).color);
            if lit_overlay == 0
            plot(d(i).R_nE_ij,(1-d(i).sigma)*o(i).UB_intensive,'--','LineWidth',1,'Color',z(i).color);
            plot(d(i).R_nE_ij,(1-d(i).sigma)*o(i).LB_intensive,'--','LineWidth',1,'Color',z(i).color);
            end
        end

        if lit_overlay == 1
            load('literature.mat')
%             a(num+1)=plot(ln_ij,(1-d(i).sigma)*(1-elast_r_pareto(range1)./elast_e_pareto(range1)),'-','LineWidth',2.5,'Color','Red');
            a(num+1)=plot(ln_ij,(1-d(i).sigma)*(1-elast_r_tpareto(range1)./elast_e_tpareto(range1)),'-','LineWidth',2.5,'Color','Green');
            a(num+2)=plot(ln_ij,(1-d(i).sigma)*(1-elast_r_logn(range1)./elast_e_logn(range1)),'-','LineWidth',2.5,'Color','Blue');
            %p4=plot(log(unif_EKK(range2)),(1-d(i).sigma)*(1-smoothdata(elast_r_EKK(range2),'gaussian',10)./smoothdata(elast_e_EKK(range2),'gaussian',10)),'-','LineWidth',2.5,'Color','Orange');
        end    

        ylabel('Intensive Margin Elasticity','FontSize',18)

    case 'extensive'
        clf
        hold on
        for i= 1:num
            a(i) = plot(d(i).R_nE_ij,(1-d(i).sigma)*1./o(i).epsilon  ,'-','LineWidth',3,'Color',z(i).color);
%             if lit_overlay == 0
            plot(d(i).R_nE_ij,-(1-d(i).sigma)*o(i).UB_extensive,'--','LineWidth',1,'Color',z(i).color);
            plot(d(i).R_nE_ij,-(1-d(i).sigma)*o(i).LB_extensive,'--','LineWidth',1,'Color',z(i).color);
%             end
        end

        if lit_overlay == 1
            load('literature.mat')
%             a(num+1)=plot(ln_ij,(1-d(i).sigma)*1./elast_e_pareto(range1),'-','LineWidth',2,'Color',[0 0.4470 0.7410]);
            a(num+1)=plot(ln_ij,(1-d(i).sigma)*1./elast_e_tpareto(range1),'-','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
            a(num+2)=plot(ln_ij,(1-d(i).sigma)*1./elast_e_logn(range1),'-','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
            a(num+3)=plot(log(unif_EKK(range2)),(1-d(i).sigma)*1./smoothdata(elast_e_EKK(range2),'gaussian',10),'-','LineWidth',2,'Color',[0.4940 0.1840 0.5560]	);
        end    


%         ylabel('Extensive Margin Elasticity','FontSize',18)
        ylabel('$\theta^{e}(\sigma-1)$','FontSize',24,'Interpreter','Latex')

    end




    hold off  
    xlabel('Log Exporter Firm Share','FontSize',18)
    ax = gca;
    ax.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];
    ax.XTickLabel = {'0.01%' ,'0.1%','1%','10%','100%'};
    ax.FontSize = 18; 
    title(gtitle);

    if lit_overlay == 1
        titles = {o.title};
%         titles{num+1} = 'Pareto';
        titles{num+1} = 'Truncated Pareto';
        titles{num+2} = 'Lognormal';
        titles{num+3} = 'EKK';
        xlim([-8 -1])
        xlim([-8 -2]);
    else
        xlim([-9.9 -2]);
    end

    switch filename
        case  'output/F2ab_dual_'
            legend('off')
        case  'output/F4b_HS_origin'

        otherwise
            if lit_overlay == 1
            else
                titles = {o.title};
                for i= 1:num
                    disp(i)
                    legend(a,titles,'FontSize',18,'Location', 'Best')
                end
            end
            disp(i)
            legend(a,titles,'FontSize',18,'Location', 'Best')
    end

    box on
    legend('boxoff')
    saveas(a(i),[filename,type],'epsc')
end

