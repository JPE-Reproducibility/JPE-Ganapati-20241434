%% Wrapper to Plot Aggregate Data
    
    a = scatter(ln_n_ij(group1),metric(group1),[],[0.6350 0.0780 0.1840],'filled','o');
    hold on
    b = scatter(ln_n_ij(group2),metric(group2),[],[0 0.4470 0.7410],'filled','square');
    text(ln_n_ij(all),metric(all),labels.textdata(all),'VerticalAlignment','top','HorizontalAlignment','left')
    hold off
    xlabel('Average Exporter Share','FontSize',16,'Interpreter','Latex')
    f =gca;
    if scenario == 2
        ylim([-.5 .35]);
    else
        ylim([-1 1.2]);
    end
    f.XTick   = -9:-3;f.XTickLabel = {'0.01%' ,'0.03%','0.09%','0.25%','0.67%','1.8%','5%'};
