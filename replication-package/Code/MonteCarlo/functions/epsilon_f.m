function e = epsilon_f(nvec, e_bar, alpha_e, gamma_e)
    
    e = zeros(size(nvec));
    for k = 1:numel(nvec)
        e(k) = generate_eps(1 - nvec(k), e_bar, alpha_e, gamma_e);
    end
end