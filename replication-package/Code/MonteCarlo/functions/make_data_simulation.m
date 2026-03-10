function [ d ] = make_data_simulation( M,knots,type)

%% Make Data File from Matrix M
    filter = find(ones(size(M(:,1)))==1);
    d.N           = size(filter,1);

%     filter        = -1>(M(:,3));
%     filter        = -1>(M(:,3)) & -1.5<(M(:,3));
%     filter        = -1>(M(:,3)) & -2<(M(:,3));
%     d.N           = sum(filter,1);

    d.I           = M(filter,1);
    d.J           = M(filter,2);
    d.R_xbar_ij   = M(filter,4);            % lx_bar
    d.R_nE_ij     = (M(filter,3));         % ln_ij
    d.ldist_raw   = (M(filter,5));             
    d.ltariff_raw = (M(filter,6));            
    d.X_ij        = (M(filter,7));            
    d.R_I         = dummyvar(M(filter,1));
    d.R_J         = dummyvar(grp2idx(M(filter,2)));

%% Create Fixed Effects Dummies
    % FE are different on J
    d.FEA = [ d.R_I d.R_I*0 d.R_J(:,1:(size(d.R_J,2)-1))   d.R_J(:,1:(size(d.R_J,2)-1))*0];
    d.FEB = [ d.R_I*0 d.R_I d.R_J(:,1:(size(d.R_J,2)-1))*0 d.R_J(:,1:(size(d.R_J,2)-1))  ];

%     d.FEA = [ d.R_J(:,1:(size(d.R_J,2)))   d.R_J(:,1:(size(d.R_J,2)))*0];
%     d.FEB = [ d.R_J(:,1:(size(d.R_J,2)))*0 d.R_J(:,1:(size(d.R_J,2)))  ];
    
%     d.FEA = [ones(size(d.I)) ones(size(d.I))*0];
%     d.FEB = [ones(size(d.I))*0 ones(size(d.I))];

%% MAKE KNOTS

    if knots == 1
        Xknots  = d.R_nE_ij;   
        d.Deriv = ones(size(Xknots));
        k       = [];
        d.kG{1} = [];
    else
        x = d.R_nE_ij;
        [Xknots,k,Deriv] = evknots(knots,x,type) ; 
        d.Deriv     = Deriv;
        d.kG{1} = k;
    end


%% Create Insturments
    d.Z1     =  d.ldist_raw;
    d.Z2     =  d.ltariff_raw;

    [Zknots1b] = evknots(knots,d.Z1,type) ; 
%     [Zknots2b] = evknots(knots,d.Z2,type) ; 

%     Zknots = [ Zknots1b Zknots2b];
    Zknots =  Zknots1b ;

%%
    d.Z = Zknots;
    d.ZA_K  = [Zknots  d.FEA  ];  
    d.ZB_K  = [Zknots  d.FEB  ];  

    % Below is some slightly sloppy coding - different subroutines use
    % different variable names
    d.Xknots    = Xknots;
    d.k         = k;
    d.X_K       = Xknots;
    d.FE_C      = [];

%% Pre-Compute Projection Matrices (to save on computing cost down the road)
    disp('create projection')
    d.Z_K       = [d.ZA_K;d.ZB_K];
    d.PZ_K      = d.Z_K*(d.Z_K'*d.Z_K)^-1*d.Z_K';
    d.FE        = [[d.FEA d.FE_C];[d.FEB d.FE_C]];
    d.FE_PZ_Ki  = (d.FE'*d.PZ_K*d.FE)^-1*d.FE'*d.PZ_K;

%     CHECK IF THIS SEEMS REASONABLE (SHOULDN'T BE TOO BIG)
    sum(sum(inv(d.Z_K'*d.Z_K)))


end


