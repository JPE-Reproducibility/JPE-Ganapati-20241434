function plot_logcorr_pareto(alpha, gamma)
% PLOTS: theta^e(n) vs ln n, CDF, and PDF for the log-corrected Pareto.
%   Survival:  S(x) = (L^alpha)(log(L)^gamma) * x^(-alpha) * (log(x))^(-gamma)
%   Hazard:    r(x) = alpha/x + gamma/(x*log(x))
%   with L = e (so S(L)=1 and F(L)=0).  x0 = 1 by construction.
%
% Inputs:
%   alpha > 0, gamma \in R
%
% Notes:
% - If alpha + gamma <= 0, the hazard at x=L may be <= 0.  In that case,
%   we lift L slightly above exp(-gamma/alpha) to keep r(x)>0 on [L, inf).
% - theta^e(n) = epsilon(n) * r(epsilon(n)), where epsilon solves S(epsilon)=n.

    arguments
        alpha (1,1) double {mustBePositive}
        gamma (1,1) double
    end

    % ----- Left endpoint and normalization -----
    L = exp(1);                       % default left endpoint
    if alpha + gamma <= 0
        % ensure r(L)>0
        L = max(L, exp(-gamma/alpha) * 1.01);
        warning('Adjusted left endpoint to L=%.6g to keep hazard positive.', L);
    end
    S  = @(x) (L.^alpha) .* (log(L).^gamma) .* x.^(-alpha) .* (log(x)).^(-gamma);
    r  = @(x) alpha./x + gamma./(x.*log(x));
    f  = @(x) r(x).*S(x);
    F  = @(x) 1 - S(x);

    % ----- θ^e(n) on n ∈ (1e-4, 1) -----
    n  = linspace(1e-4, 0.9999, 1500);
    lnn = log(n);
    eps_n = arrayfun(@(p) inv_tail_quantile(p, alpha, gamma), n);
    theta = eps_n .* r(eps_n);               % θ^e(n) = ε r(ε)

%     figure('Color','w'); 
%     plot(lnn, theta, 'LineWidth', 1.8);
%     xlabel('ln n'); ylabel('\theta^e(n)');
%     title(sprintf('\\theta^e(n) vs ln n  (\\alpha=%.3g, \\gamma=%.3g, L=%.3g)', alpha, gamma, L));
%     grid on;

    % ----- PDF: wide log–log view -----
    x1 = linspace(L, 200, 2000);
%     figure('Color','w');
%     loglog(x1, max(f(x1), 0), 'LineWidth', 1.8);
%     xlabel('x'); ylabel('PDF');
%     title(sprintf('PDF (log–log),  \\alpha=%.3g, \\gamma=%.3g, L=%.3g', alpha, gamma, L));
%     grid on;

    % ----- CDF & body PDF on [L, 10] -----
    x2max = max(10, L*1.01);                   % in case L>10
    x2 = linspace(L, x2max, 1200);

%     figure('Color','w');
%     plot(x2, F(x2), 'LineWidth', 1.8);
%     xlabel('x'); ylabel('CDF');
%     title(sprintf('CDF on [L, 10]  (\\alpha=%.3g, \\gamma=%.3g, L=%.3g)', alpha, gamma, L));
%     grid on;

    figure('Color','w');
    plot(x2, max(f(x2), 0), 'LineWidth', 1.8);
    xlabel('x'); ylabel('PDF');
    title(sprintf('PDF on [L, 10]  (\\alpha=%.3g, \\gamma=%.3g, L=%.3g)', alpha, gamma, L));
    grid on;
    ax = gca;
    ax.FontSize = 16; 

    % ===== nested: tail-quantile ε(n) via solve in t = ln ε =====
    function eps = inv_tail_quantile(p, a, g)
        % Solve S(e^t) = p  =>  a - a t - g ln t = ln p
        % Monotone in t, unique solution with t > ln L >= 1.
        ln_p = log(p);
        obj  = @(t) a - a*t - g*log(t) - ln_p;
        tL   = max(1 + 1e-10, log(L));      % lower bracket (t >= ln L)
        % heuristic upper bracket that grows with |ln p|
        tU   = max(tL + 1, 1 - ln_p/a + 10);
        % fzero with bracket [tL, tU]
        t = fzero(obj, [tL, tU]);
        eps = exp(t);
    end
end
