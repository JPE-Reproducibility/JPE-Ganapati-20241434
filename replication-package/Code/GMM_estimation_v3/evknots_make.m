function [Xknots,k,Deriv]=evknots_make(k,x) 
% Make a matri to apply evknots

    V = zeros(length(x),length(k)-1);
    Deriv = V;

    V(:,1)=x;
    Deriv(:,1)=1;

    n=length(k);
    POS = @(x) x.*(x>0);
    
    for i=1:(n-2)
        NUM1 = POS((x-k(n-1))).^3 *(k(n)  -k(i));
        NUM2 = POS((x-k(n)  )).^3 *(k(n-1)-k(i));
        NUM  = POS((x-k(i)  )).^3 -(k(n)-k(n-1))^(-1)*(NUM1-NUM2);
        DEN = (k(n)-k(1))^2;
        V(:,(i+1))=NUM./DEN;
    end
    for i=1:(n-2)
        NUM1 = 3*POS((x-k(n-1))).^2 *(k(n)  -k(i));
        NUM2 = 3*POS((x-k(n)  )).^2 *(k(n-1)-k(i));
        NUM  = 3*POS((x-k(i)  )).^2 -(k(n)-k(n-1))^(-1)*(NUM1-NUM2);
        DEN = (k(n)-k(1))^2;
        Deriv(:,(i+1))=NUM./DEN;
    end
    Xknots = V;
    
end

