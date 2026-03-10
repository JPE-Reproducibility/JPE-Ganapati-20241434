function [Xknots,kG,Deriv]=evknots_makeG(kG,x,G,gg2) 
% Make a matri to apply evknots

    V = zeros(length(x),length(kG{1})-1);
    Deriv = V;

    if size(x,2) > 1
        x = x';
        G = G';
    end

    V(:,1)=x;
    Deriv(:,1)=1;

    if size(x,1) > size(G,1)
        G = repmat(G,size(x));
    end

    for g = 1:size(kG,2)

        n=length(kG{g});
        POS = @(x) x.*(x>0);
        k = kG{g};

        gg = g-1;
%         display(gg)
        for i=1:(n-2)
            NUM1 = POS((x(G==gg)-k(n-1))).^3 *(k(n)  -k(i));
            NUM2 = POS((x(G==gg)-k(n)  )).^3 *(k(n-1)-k(i));
            NUM  = POS((x(G==gg)-k(i)  )).^3 -(k(n)-k(n-1))^(-1)*(NUM1-NUM2);
            DEN = (k(n)-k(1))^2;
            V((G==gg),(i+1))=NUM./DEN;
        end
        for i=1:(n-2)
            NUM1 = 3*POS((x(G==gg)-k(n-1))).^2 *(k(n)  -k(i));
            NUM2 = 3*POS((x(G==gg)-k(n)  )).^2 *(k(n-1)-k(i));
            NUM  = 3*POS((x(G==gg)-k(i)  )).^2 -(k(n)-k(n-1))^(-1)*(NUM1-NUM2);
            DEN = (k(n)-k(1))^2;
            Deriv((G==gg),(i+1))=NUM./DEN;
        end
    end
    Xknots = V;

%     if nargin < 4
%     %     Convert to index
%         dc0 = repmat(G==0,1,max(1,length(kG{1})-1));
%         dc1 = repmat(G==1,1,max(1,length(kG{1})-1));
%         dc2 = repmat(G==2,1,max(1,length(kG{1})-1));
%         dc3 = repmat(G==3,1,max(1,length(kG{1})-1));
%     
%         Xknots =   [ Xknots.*dc0  Xknots.*dc1  Xknots.*dc2  Xknots.*dc3];
%         Deriv = [  Deriv.*dc0   Deriv.*dc1   Deriv.*dc2   Deriv.*dc3];
%     else
%         dc0 = repmat(G==0,1,max(1,length(kG{1})-1));
%         dc1 = repmat(G==1,1,max(1,length(kG{1})-1));
%     
%         Xknots =   [ Xknots.*dc0  Xknots.*dc1 ];
%         Deriv = [  Deriv.*dc0   Deriv.*dc1   ];

        Xknots2 = [];
        Deriv2 = [];
        for i=1:max(size(kG))
            dc = repmat(G==(i-1),1,max(1,length(kG{i})-1));
            Xknots2 = [Xknots2 Xknots.*dc];
            Deriv2 = [ Deriv2 Deriv.*dc];
        end
        Xknots = Xknots2;
        Deriv = Deriv2;
%     end



end

