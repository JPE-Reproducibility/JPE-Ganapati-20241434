function [Xknots,k,Deriv]=evknots(knots_points,x,type) 
%     knots=prctile(x,linspace(0,100,knots_points))
    if knots_points == 2
        k=prctile(x,[10  90 ]) % From Harrell (2001, 23)

    elseif knots_points == 3
        k=prctile(x,[10 50 90 ]) % From Harrell (2001, 23)
    elseif knots_points == 33
        k=prctile(x,[25  50 75 ]) 
    elseif knots_points == 34
        k=prctile(x,[25  50 70 ]) 
    elseif knots_points == 4
        k=prctile(x,[5 35 65 95]) % From Harrell (2001, 23)
    elseif knots_points == 5
        k=prctile(x,[5 27.5 50 72.5 95]) % From Harrell (2001, 23)
            elseif knots_points == 6
        k=prctile(x,[5 23 41 59 77 95]) % From Harrell (2001, 23)
            elseif knots_points == 7
        k=prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5]) % From Harrell (2001, 23)
    else
        k=prctile(x,linspace(0,100,knots_points+2));
        k=k(2:(knots_points+1))
    end
    V = zeros(length(x),length(k)-1);
    Deriv = V;
    
    POS = @(x) x.*(x>0);
    
    switch lower(type)
        
        case 'linear'
%             NOT DEBUGGED - DO NOT USE YET
            V = [];
            Deriv = [];
            n=length(k)+1;
            
            V(:,1)=x.*(x<=k(1));
            for i = 2:(n-1)
                V(:,i) = x.*(x>k(i-1)).*(x<=k(i));
            end
            V(:,n) = x.*(x>k(n-1));
            
            Deriv(:,1) = (x<=k(1));
            for i = 2:(n-1)
                Deriv(:,i) = (x>k(i-1)) & (x<=k(i));
            end
            Deriv(:,n) = (x>k(n-1));
            Xknots = V;          
        case 'lineara'
%             NOT DEBUGGED - DO NOT USE YET
            V = [];
            Deriv = [];
            n=length(k)+1;

            V(:,1)=min(x,k(1));
            for i = 2:(n-1)
                V(:,i) = max(min(x,k(i)),k(i-1)) - k(i-1);
            end
            V(:,n) = max(x,k(n-1))-k(n-1);
            
            Deriv(:,1) = (x<=k(1));
            for i = 2:(n-1)
                Deriv(:,i) = (x>k(i-1)) & (x<=k(i));
            end
            Deriv(:,n) = (x>k(n-1));
            Xknots = V;
        otherwise

            V(:,1)=x;
            Deriv(:,1)=1;

            n=length(k);
    
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


end
