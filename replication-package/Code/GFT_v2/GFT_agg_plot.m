%% Wrapper to Plot Aggregate Data
    all = poor | rich;
    a = scatter(ln_n_ij(rich),metric(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    b = scatter(ln_n_ij(poor),metric(poor),[],[0 0.4470 0.7410],'filled','square');
    text(ln_n_ij(all),metric(all),labels(all),'VerticalAlignment','top','HorizontalAlignment','right')
    hold off
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
    
    f =gca;
    ylim([-.4 .6]);
    legend([a b ],{'Developed', 'Developing'},'FontSize',14,'Location','best')
    f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
