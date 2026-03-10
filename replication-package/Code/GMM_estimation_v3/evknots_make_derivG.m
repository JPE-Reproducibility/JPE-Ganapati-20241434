function [Deriv]=evknots_make_derivG(kG,x,G,gg2)
% Make a matri to apply evknots

V = zeros(length(x),length(kG{1})-1);
V(:,1)=1;

if size(x,2) > 1
    x = x';
    G = G';
end

if size(x,1) > size(G,1)
    G = repmat(G,size(x));
end

for g = 1:size(kG,2)

    n=length(kG{g});
    POS = @(x) x.*(x>0);
    k = kG{g};
    gg = g-1;
    for i=1:(n-2)
        NUM1 = 3*POS((x(G==gg)-k(n-1))).^2 *(k(n)  -k(i));
        NUM2 = 3*POS((x(G==gg)-k(n)  )).^2 *(k(n-1)-k(i));
        NUM  = 3*POS((x(G==gg)-k(i)  )).^2 -(k(n)-k(n-1))^(-1)*(NUM1-NUM2);
        DEN = (k(n)-k(1))^2;
        V((G==gg),(i+1))=NUM./DEN;
    end

    Deriv = V;

end

Deriv2 = [];
for i=1:max(size(kG))
    dc = repmat(G==(i-1),1,max(1,length(kG{i})-1));
    Deriv2 = [ Deriv2 Deriv.*dc];
end
Deriv = Deriv2;

end

