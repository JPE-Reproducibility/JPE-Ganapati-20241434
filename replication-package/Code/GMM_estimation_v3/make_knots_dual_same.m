d.kG{1} = []; d.kG{2} = [];

if knots>2
    x = d.R_nE_ij;

    [Xknots1t,k1t,Deriv1t] = evknots(knots,x,type) ;
    k1a = k1t;
    k1b = k1t;
    Xknots1a =Xknots1t(d.category==0,:);
    Xknots1b =Xknots1t(d.category==1,:);
    Deriv1a =Deriv1t(d.category==0,:);
    Deriv1b =Deriv1t(d.category==1,:);

%     [Xknots1a,k1a,Deriv1a] = evknots(knots,x(d.category==0),type) ;
%     [Xknots1b,k1b,Deriv1b] = evknots(knots,x(d.category==1),type) ;

    d.k = [k1a k1b  ];
    d.k1a = k1a;    d.k1b = k1b;
    d.kG{1} = k1a;    d.kG{2} = k1b;

    Xknots(d.category==0,:)=Xknots1a;
    Xknots(d.category==1,:)=Xknots1b;

    Deriv(d.category==0,:)=Deriv1a;
    Deriv(d.category==1,:)=Deriv1b;
    d.Deriv     = Deriv;

    z1a       = d.ldist_raw(d.category==0);
    z1b       = d.ldist_raw(d.category==1);
    
    ztemp1a   = zeros(length(Xknots1a), 1);
    ztemp1b   = zeros(length(Xknots1b), 1);

    for i   = 1:length(z1a)
        ztemp1a(i) = sum(z1a <= z1a(i))/length(z1a);
    end
    for i   = 1:length(z1b)
        ztemp1b(i) = sum(z1b <= z1b(i))/length(z1b);
    end

    [Zknots1a, ~] = cosine_basis(ztemp1a, 11);
    [Zknots1b, ~] = cosine_basis(ztemp1b, 11);

    z2a       = d.ltariff_raw(d.category==0);
    z2b       = d.ltariff_raw(d.category==1);
    ztemp2a   = zeros(length(Xknots1a), 1);
    ztemp2b   = zeros(length(Xknots1b), 1);
    for i   = 1:length(z2a)
        ztemp2a(i) = sum(z2a <= z2a(i))/length(z2a);
    end
    for i   = 1:length(z2b)
        ztemp2b(i) = sum(z2b <= z2b(i))/length(z2b);
    end
    [Zknots2a, ~] = cosine_basis(ztemp2a, 11);
    [Zknots2b, ~] = cosine_basis(ztemp2b, 11);

    t = [Zknots1a(:,2:end) Zknots2a(:,2:end)];
    Zknots = zeros(size(Zknots,1),size(t,2));
    Zknots(d.category==0,:)=[Zknots1a(:,2:end) Zknots2a(:,2:end)];
    Zknots(d.category==1,:)=[Zknots1b(:,2:end) Zknots2b(:,2:end)];
end