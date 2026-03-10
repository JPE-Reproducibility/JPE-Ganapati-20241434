function [F,elastlnbeta] = match_LogPareto(e_high,e_low,alpha,beta,gamma,elast_epsd,lnn)

    e_beta = @(n) logpareto_inv_shifted(1-n,alpha,beta,gamma);

    lne_beta = @(lnn) log(e_beta(exp(lnn)));

    dlne = diff(lne_beta(lnn));

    DD=mean(diff(lnn));

    dlne(lnn<e_low)  = min(elast_epsd)*DD;
    dlne(lnn>e_high) = max(elast_epsd)*DD;
    N = max(size(elast_epsd));

    elastlnbeta = dlne(1:N)./diff(lnn);
    [~, index_low] = min(abs(lnn - e_low));
    [~, index_high] = min(abs(lnn - e_high));
    A = 10000*abs(elastlnbeta((index_low))-min(elast_epsd))^2;
    B = 10000*abs(elastlnbeta((index_high))-max(elast_epsd))^2;
    F = sum((elast_epsd(1:N)-elastlnbeta).^2) + A +B;

    if alpha<0 || beta <0
        F = F + 500000*alpha^2 + 500000*beta^2;
    else
    end
    
end