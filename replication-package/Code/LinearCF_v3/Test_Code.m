
%% Testing
if testing == 111

    inner = zeros(dim);
    for i = 1:dim
        for j = 1:dim
            for d = 1:dim
                inner(i,j) = inner(i,j) + y_ij_init(i,d).*x_ij_init(j,d);
            end
        end
    end
    
    
    dlnr = log(r_2./r_1);
    
    left = zeros(1,dim)
    for j = 1:dim
        for o = 1:dim
                left(1,j) = left(1,j) + x_ij_init(o,j).*Gtheta(0,0).*dlnr(o,j);
        end
    end
    
    w_guess = zeros(dim,1);
    
    FFF = @(w_guess) theta(0)*(sigma).*w_guess-sum((y_ij_init  + (Gtheta(0,0)*sigma- 1).*inner).*w_guess',2) ...
        -sum(y_ij_init .* (Gtheta(0,0).*dlnr - left),2);
    
    start = [w_guess];
    solved = fsolve(@(ww) FFF([ww]),0*ones(dim,1),optimoptions('fsolve','Display','iter','FunctionTolerance',1e-12));
    dlw_simple = [solved];
    
    % with normalization
    solved = fsolve(@(ww) FFF([ww;0]),0*ones(dim-1,1),optimoptions('fsolve','Display','iter'));
    dlw_simple = [solved;0];
    
    FFF(dlw_simple);
    % FFF([0 ; 0])
    FFF(log(dw_linear));
    
    % dlw_simple = [3 ; 3]
    
    dlP_simple = ((1-theta(0)).*dlw_simple' ...
        -sum(x_ij_init.*(theta(0).*dlnr+(1-Gtheta(0,0)*sigma).*dlw_simple),1)) ...
        ./(theta(0)*(sigma-1));
    
    dW_simple   = exp(dlw_simple)./exp(dlP_simple');
    
    
    [log(dw_linear) dlw_simple ];
    [sum(dlnP_sequence,2) dlP_simple'];
    [log(Real_Wage) log(dW_simple)];
    
    
    % Compares Real wages
    % log(Real_Wage)./log(dW_simple)
    % scatter(log(Real_Wage), log(dRW))
    % scatter(log(dw_linear), dlw_linear_simple)
    
    
%% Loops if Pareto with one segment
    if segments == 1

        dlX = theta(0).*dlnr  ...
            + (1-theta(0)*sigma)*dlw_simple ...
            + theta(0)*(sigma-1)*(dlP_simple) ...
            + theta(0) .* dlw_simple';
                
        dlN_simple = sum(y_ij_init.*(theta(0)-1)*(-dlnr+sigma*dlw_simple-(sigma-1)*dlP_simple-dlw_simple'),2);
        

        dlN_t = sum(y_ij_init.*(theta(0)-1)*(-dlnr+sigma*sum(dlnw_sequence,2)-(sigma-1)*sum(dlnP_sequence,2)'-sum(dlnw_sequence,2)'),2);
        
        dln_simple = 1/elast_epsilon(0).*(-dlnr+sigma*dlw_simple-(sigma-1)*dlP_simple);
        
        dlnE =  (dlnw );
        
        dlnn = (-dln_r ...
            + repmat(sigma*dlnw,1,dim) ...
            - repmat((sigma-1)*log(dP_linear)',dim,1) ...
            - repmat(dlnE',dim,1)) ...
        ./Gelast_epsilon(n_ij,G_ij)    ;
        sum(x_ij_init,1);
        
        
        dlnN = sum(y_ij_init.*theta(0).*elast_r(0).*dlnn./(sigma-1),2);
        dlnN_simple2 = sum(y_ij_init.*theta(0).*elast_r(0).*dln_simple./(sigma-1),2);
        
        
        dlXX = dlnn_sequence+dlnxbar_sequence;
        
        sum(x_ij_init.* dlXX,1);
        sum(y_ij_init.* dlXX,2);
        
        
        sum(y_ij_init.*dlnn_sequence,2);
        sum(y_ij_init.*dln_simple,2);
    end
    
end




