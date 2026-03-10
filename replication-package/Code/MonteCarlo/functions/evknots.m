function [Xknots,k,Deriv]=evknots(knots_points,x,type) 
    if knots_points == 2
        k=prctile(x,[10  90 ]); % From Harrell (2001, 23)
    elseif knots_points == 3
        k=prctile(x,[10 50 90 ]); % From Harrell (2001, 23)
    elseif knots_points == 33
        k=prctile(x,[25  50 75 ]); 
    elseif knots_points == 34
        k=prctile(x,[25  50 70 ]); 
    elseif knots_points == 4
        k=prctile(x,[5 35 65 95]); % From Harrell (2001, 23)
    elseif knots_points == 5
        k=prctile(x,[5 27.5 50 72.5 95]); % From Harrell (2001, 23)
    elseif knots_points == 55
    %  FOR LOGNORMAL

%         k=prctile(x,[5 27.5 50  72.5 95 98 98.25 98.5 99 99.5]); % Works
%         k=prctile(x,[2.5  65.83 98 98.25 98.5    99.75]); % A bit iffy at
%         the edges
        % k=prctile(x,[2.5 18.33 34.17 50 65.83 98 98.25 98.5    99.75]); % Works
                % % k=prctile(x,[1 5 27.5 50 72.5 95 99.5 99.75]); % From Harrell (2001, 23)
        k=[prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 99])]; % Old Draft?
        % k=prctile(x,[10 20 30 10 50 60 70 80 90 95 98 99 100]); % From Harrell (2001, 23)
        % k=[prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 99])]; % From Harrell (2001, 23)
        % k=[prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 99 ])]; % From Harrell (2001, 23)
        k=[prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 ])]; % From Harrell (2001, 23)



        % Nov 4 play
        % k=[prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 ])]; % From Harrell (2001, 23)
        % k=[prctile(x,[2.5 18.33 34.17 50  81.67 97.5 ])]; % From Harrell (2001, 23)
        k=[prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 99])]; % Old Draft?
        k=prctile(x,[2.5 5 27.5 50 72.5 97.5 99]); % From Harrell (2001, 23)
        k=prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5         98 99 99.5 ]); % I made Up
        k=prctile(x,[2.5 18.33 34.17 50 81.67 97.5         98 99 99.5 99.9]); % I made Up
        k=prctile(x,[2.5 18.33 34.17 50 97.5         98 99 99.5 99.9]); % I made Up
        exp(k)*100

    elseif knots_points == 56
%         k=prctile(x,[2.5 18.33 34.17 50 65.83 98 98.25 98.5    99.75]); % Works
%         k=prctile(x,[5 27.5 50 72.5 97 98 99  ]); % From Harrell (2001, 23)
%         k=prctile(min(x+rand(size(x))*.002,max(x)),[2.5 18.33 34.17 50 65.83 98 98.25 98.4 98.5  98.6  99.75]) % Works
%         k=prctile(x,[2.5 18.33 34.17 50 65.83 98 98.25 98.5    99.75]); % Works
%         k=prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 98.5 99.8]); % From Harrell (2001, 23)
%         k=prctile(x,[2.5 18.33 34.17 50 65.83 81.67 85 90 95 96 97.5 98.5 99.8]); % From Harrell (2001, 23)

        k=prctile(x,[2.5 18.33 34.17 50 65.83 90 98 98.25   99 99.75]); % Works
        k=prctile(x,[5 27.5 50 72.5 99.6 99.75 99.9 99.95]); % From Harrell (2001, 23)
        k=prctile(x,[5 27.5 50 72.5  99.75 99.9 99.95]); % From Harrell (2001, 23)


        % Nov 4 play
        k=[prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5 99])]; % Old Draft?


    elseif knots_points == 57
        % k=prctile(x,[5 27.5  66 95]) % From Harrell (2001, 23)
        % k=[-9.88055994 -8	-7.210481 -5	-3.66798696 ]; % hardcoded from the paper
        % k=[-9.88055994 -8.5	-7.810481 -5	-3.66798696 ]; % hardcoded from the paper
        % [-9.8806 -7.0210 -3.6680]
        % k=prctile(x,[5 35 76 78 95]); % From Harrell (2001, 23)
        % k=prctile(x,[5 20 63 66 95]) % From Harrell (2001, 23)
        % k=prctile(x,[5 20 63 66 89]) % Ad hoc; built from support
        % k=prctile(x,[5 20 63 67 89]) % Ad hoc; built from support
        k=prctile(x,[5 20 65 67 89]) % Ad hoc; built from support
        % k=prctile(x,[5 20 65 67 89]+1) % Ad hoc; built from support
        % k=prctile(x,[5 15 63 67 89]) % Ad hoc; built from support
        % k=prctile(x,[5 20 65 67 90]) % Ad hoc; built from support
    elseif knots_points == 58
        k=prctile(x,[0.5 18.33 34.17 50 97.5         98 99 99.5 99.9]); % ST nov 5 play
        

    elseif knots_points == 6
        k=prctile(x,[5 23 41 59 77 95]); % From Harrell (2001, 23)
    elseif knots_points == 7
        k=prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5]); % From Harrell (2001, 23)
    elseif knots_points == 77
%         k=prctile(x,[2.5 18.33 34.17 50 65.83 81.67 90 95 97.5 98 99 99.5 ]); % Made up, works for GFT
        k=prctile(x,[2.5 18.33 34.17 50 65.83 81.67 97.5         98 99 99.5 ]); % I made Up
    elseif knots_points == 9
        k=prctile(x,[1 2.5 18.33 34.17 50 65.83 81.67 97.5 99 99.25 99.5 99.7 99.8]); % I made Up
    elseif knots_points == 99
        k=[-9.88055994	-7.0210481	-3.66798696 ]; % hardcoded from the paper
    else
        k=prctile(x,linspace(0,100,knots_points+2));
        k=k(2:(knots_points+1));
    end
    V = zeros(length(x),length(k)-1);
    Deriv = V;
    
    POS = @(x) x.*(x>0);
    
    switch lower(type)
        
        case 'linear'
%             NOT DEBUGGED - DO NOT USE YET
        case 'lineara'
%             NOT DEBUGGED - DO NOT USE YET
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
