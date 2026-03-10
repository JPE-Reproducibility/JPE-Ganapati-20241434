function ge = generate_H(e, e_bar, alpha_e, gamma_e)
    % generates 
    % G^e(e) = 1 - (e/e_bar)^(-alpha_e) * (ln(e)/ln(e_bar))^(-gamma_e)
    % for a given e

    ge = 1 - (e./e_bar).^(-alpha_e) .* (log(e)./log(e_bar)).^(-gamma_e);
end