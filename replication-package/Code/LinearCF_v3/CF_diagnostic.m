%% Wrapper to Plot Diagnostic Data
    a = scatter(ln_n_ij(rich),plot_var(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    b = scatter(ln_n_ij(poor),plot_var(poor),[],[0 0.4470 0.7410],'filled','square');
    text(ln_n_ij(all),plot_var(all),labels.textdata(all),'VerticalAlignment','top','HorizontalAlignment','right')
    hold off
    f =gca;
    f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
    % xlabel('Average Exporter Share $\log {\bar{n}_{ij}}$','FontSize',16,'Interpreter','Latex')
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
    legend([a b ],{'Developed Countries', 'Developing Countries'},'FontSize',14,'Location','best')