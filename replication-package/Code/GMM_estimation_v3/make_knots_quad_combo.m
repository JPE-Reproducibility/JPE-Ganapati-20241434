d.kG{1} = []; d.kG{2} = [];
d.kG{3} = []; d.kG{4} = [];

if knots>2
    x = d.R_nE_ij;
    [Xknots1a,k1a,Deriv1a] = evknots(knots,x(d.category==0),type) ;
    [Xknots1b,k1b,Deriv1b] = evknots(knots,x(d.category==1),type) ;
    [Xknots1c,k1c,Deriv1c] = evknots(knots,x(d.category==2),type) ;
    [Xknots1d,k1d,Deriv1d] = evknots(knots,x(d.category==3),type) ;
    d.k = [k1a k1b k1c k1d];

    d.k1a = k1a;    d.k1b = k1b;
    d.k1c = k1c;    d.k1d = k1d;

    d.kG{1} = k1a;    d.kG{2} = k1b;
    d.kG{3} = k1c;    d.kG{4} = k1d;

    Xknots(d.category==0,:)=Xknots1a;
    Xknots(d.category==1,:)=Xknots1b;
    Xknots(d.category==2,:)=Xknots1c;
    Xknots(d.category==3,:)=Xknots1d;

    Deriv(d.category==0,:)=Deriv1a;
    Deriv(d.category==1,:)=Deriv1b;
    Deriv(d.category==2,:)=Deriv1c;
    Deriv(d.category==3,:)=Deriv1d;

    d.Deriv     = Deriv;
    d.Deriv_dim = Deriv;

    z1a       = d.ldist_raw(d.category==0);
    z1b       = d.ldist_raw(d.category==1);
    z1c       = d.ldist_raw(d.category==2);
    z1d       = d.ldist_raw(d.category==3);

    ztemp1a   = zeros(length(Xknots1a), 1);
    ztemp1b   = zeros(length(Xknots1b), 1);
    ztemp1c   = zeros(length(Xknots1c), 1);
    ztemp1d   = zeros(length(Xknots1d), 1);

    for i   = 1:length(z1a)
        ztemp1a(i) = sum(z1a <= z1a(i))/length(z1a);
    end
    for i   = 1:length(z1b)
        ztemp1b(i) = sum(z1b <= z1b(i))/length(z1b);
    end
    for i   = 1:length(z1c)
        ztemp1c(i) = sum(z1c <= z1c(i))/length(z1c);
    end
    for i   = 1:length(z1d)
        ztemp1d(i) = sum(z1d <= z1d(i))/length(z1d);
    end

    [Zknots1a, ~] = cosine_basis(ztemp1a, 11);
    [Zknots1b, ~] = cosine_basis(ztemp1b, 11);
    [Zknots1c, ~] = cosine_basis(ztemp1c, 11);
    [Zknots1d, ~] = cosine_basis(ztemp1d, 11);

    z2a       = d.ltariff_raw(d.category==0);
    z2b       = d.ltariff_raw(d.category==1);
    z2c       = d.ltariff_raw(d.category==2);
    z2d       = d.ltariff_raw(d.category==3);

    ztemp2a   = zeros(length(Xknots1a), 1);
    ztemp2b   = zeros(length(Xknots1b), 1);
    ztemp2c   = zeros(length(Xknots1c), 1);
    ztemp2d   = zeros(length(Xknots1d), 1);

    for i   = 1:length(z2a)
        ztemp2a(i) = sum(z2a <= z2a(i))/length(z2a);
    end
    for i   = 1:length(z2b)
        ztemp2b(i) = sum(z2b <= z2b(i))/length(z2b);
    end
    for i   = 1:length(z2c)
        ztemp2c(i) = sum(z2c <= z2c(i))/length(z2c);
    end
    for i   = 1:length(z2d)
        ztemp2d(i) = sum(z2d <= z2d(i))/length(z2d);
    end

    [Zknots2a, ~] = cosine_basis(ztemp2a, 11);
    [Zknots2b, ~] = cosine_basis(ztemp2b, 11);
    [Zknots2c, ~] = cosine_basis(ztemp2c, 11);
    [Zknots2d, ~] = cosine_basis(ztemp2d, 11);

    t = [Zknots1a(:,2:end) Zknots2a(:,2:end)];
    Zknots = zeros(size(Zknots,1),size(t,2));
    Zknots(d.category==0,:)=[Zknots1a(:,2:end) Zknots2a(:,2:end)];
    Zknots(d.category==1,:)=[Zknots1b(:,2:end) Zknots2b(:,2:end)];
    Zknots(d.category==2,:)=[Zknots1c(:,2:end) Zknots2c(:,2:end)];
    Zknots(d.category==3,:)=[Zknots1d(:,2:end) Zknots2d(:,2:end)];

end