%% Wrapper to Plot Flow Data
    clf
    hold on
        scatter(n(IND==0),(plot_var(IND==0)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C1,'MarkerFaceColor',C1);   a = scatter([],[],'MarkerEdgeColor',C1,'MarkerFaceColor',C1);%a.Color = C1;%[0 0.4470 0.7410];
        scatter(n(IND==1),(plot_var(IND==1)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C1);                        b = scatter([],[],'MarkerEdgeColor',C2);%b.Color = C2;%[0.8500 0.3250 0.0980];
        scatter(n(IND==2),(plot_var(IND==2)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C3,'MarkerFaceColor',C3);   c = scatter([],[],'MarkerEdgeColor',C3,'MarkerFaceColor',C3);%c.Color = C3;%[0.9290 0.6940 0.1250];
        scatter(n(IND==3),(plot_var(IND==3)),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor',C4);                        d = scatter([],[],'MarkerEdgeColor',C4);%d.Color = C4;%[0.4940 0.1840 0.5560];
    hold off
    legend('',type0,'', type1,'', type2,'', type3,'FontSize',14,'Location','northeast')
    xlabel('$\log {n_{ij}^0}$','FontSize',16,'Interpreter','Latex')
    f =gca; f.XTick   = [-9.210340372 -6.907755279  -4.605170186  -2.302585093  0 ];f.XTickLabel = {'0.01%' ,'0.1%','1%','10%','100%'};
