function [Xknots]=evknots_make_deriv(k,x) 
% Make a matri to apply evknots
    
    V = zeros(length(x),length(k)-1);
    V(:,1)=1;
    n=length(k);
    POS = @(x) x.*(x>0);
    
    for i=1:(n-2)
        NUM1 = 3*POS((x-k(n-1))).^2 *(k(n)  -k(i));
        NUM2 = 3*POS((x-k(n)  )).^2 *(k(n-1)-k(i));
        NUM  = 3*POS((x-k(i)  )).^2 -(k(n)-k(n-1))^(-1)*(NUM1-NUM2);
        DEN = (k(n)-k(1))^2;
        V(:,(i+1))=NUM./DEN;
    end
    
    Xknots = V;
    
end

