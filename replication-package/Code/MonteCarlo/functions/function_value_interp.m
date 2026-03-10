function out = function_value_interp(x,param,grid_y)
            [xu,~,idx] = unique(x);               % xu is sorted unique(x)
            param_m = accumarray(idx, param, [], @mean);
            out = interp1(xu, param_m, grid_y, 'linear', 'extrap');
end