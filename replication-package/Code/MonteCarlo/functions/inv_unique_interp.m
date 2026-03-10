function q = inv_unique_interp(lnH, x_grid, p_grid, nvec, rbvec)
        
        xg   = x_grid(nvec, rbvec);                    % Kxpp
        Nijmat = nvec .* ones(size(xg));              % Kxpp
        RBmat   = rbvec   .* ones(size(xg));              % Kxpp
        %CDF = lnH(xg, n.*ones(size(xg)), R_bar.*ones(size(xg)));
     
        K  = numel(nvec);
        pp = numel(p_grid);
        Q  = zeros(K, pp);

        %[unique_cdf, ia, ~] = unique(CDF);
        %unique_x = xg(ia);

        
        CDFs   = lnH(xg, Nijmat, RBmat);  % ppx1
        for k = 1:K
            %fprintf('Computing: %d\n', k);
            xrow = xg(k,:).';                  % ppx1
            Fx = CDFs(k,:).';
            [Fxu, ia] = unique(Fx, 'stable');   % drop flats, keep order
            Q(k,:) = interp1(Fxu, xrow(ia), p_grid(:), 'spline', 'extrap').'; % 1xpp
        end
        q = Q;
end