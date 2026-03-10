%% Wrapper to Plot Aggregate Data
    all = poor | rich;   
    a  = scatter(xmetric(rich),ymetric1(rich),[],[0.6350 0.0780 0.1840],'filled','o');
    text(xmetric(rich),ymetric1(rich),labels(rich),'VerticalAlignment','top','HorizontalAlignment','right')
    hold on
    b  = scatter(xmetric(poor),ymetric1(poor),[],[0 0.4470 0.7410],'filled','square');
    text(xmetric(poor),ymetric1(poor),labels(poor),'VerticalAlignment','bottom','HorizontalAlignment','left')
    hold off
    legend([ a b ],{'Constant Elasticity - All Countries' 'Semiparametric - Developed', 'Semiparametric - Developing' },'FontSize',14,'Location','northeast')
    f = gcf;
