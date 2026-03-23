function [ d ] = make_data_combo( M,knots,type,interaction,subset,bs,elements)

%% Make Data File from Matrix M
    switch nargin
        case 5
            switch subset
                case 'no_n_ii'
                    filter = find(round(M(:,6))~=round(M(:,7)));
                otherwise
                    filter = find(ones(size(M(:,13)))==1);
            end
        otherwise
        filter = find(ones(size(M(:,13)))==1);
    end

    if nargin >= 6
        if bs > 0
            rng(bs,'philox')
            filter = randi([min(filter) max(filter)],max(size(filter)),1);
        end
    end


    d.M         = M;
    d.N         = size(filter,1);
    d.I         = M(filter,1);
    d.J         = M(filter,2);
    d.R_I       = dummyvar(M(filter,1));
    d.R_J       = dummyvar(grp2idx(M(filter,2)));
    d.R_N_ij    = M(filter,3);
    d.R_xbar_ij = M(filter,4);

    d.lGDPc_I   = log(M(filter,6));
    d.lGDPc_J   = log(M(filter,7));
    d.rich_I    = (M(filter,8));
    d.rich_J    = (M(filter,9));
    d.R_nE_ij   = (M(filter,10));
    d.lGDP_I    = log(M(filter,11));
    d.lGDP_J    = log(M(filter,12));
    d.colony    = (M(filter,13));
    d.contig    = (M(filter,14));
    d.fta_wto   = (M(filter,15));
    d.comcur    = (M(filter,16));
    d.comlang   = (M(filter,17));
    d.IV_dist   = (M(filter,18));
    d.IV_tarr   = (M(filter,19));
    d.lpop_od   =(M(filter,20));
    d.ldist_od  = (M(filter,21));
    d.ldist_raw = (M(filter,22));
    d.ltariff_raw = (M(filter,23));

%% Create Fixed Effects Dummies
    % FE are different on J
    d.FEA = [ d.R_I d.R_I*0 d.R_J(:,1:(size(d.R_J,2)-1))   d.R_J(:,1:(size(d.R_J,2)-1))*0];
    d.FEB = [ d.R_I*0 d.R_I d.R_J(:,1:(size(d.R_J,2)-1))*0 d.R_J(:,1:(size(d.R_J,2)-1))  ];
    
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
    [Zknots2b] = evknots(knots,d.Z2,type) ; 


    if nargin >= 6
        switch elements
            case 'dist'
                Zknots = [ Zknots1b Zknots2b];
            case 'all'
                Zknots = [ Zknots1b Zknots2b d.fta_wto d.comcur d.comlang d.colony];
            otherwise
                Zknots = [ Zknots1b Zknots2b d.fta_wto d.comcur d.comlang d.colony];
        end
    else
        Zknots = [ Zknots1b Zknots2b d.fta_wto d.comcur d.comlang d.colony];
    end



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

%% Create Interaction Effects
if nargin >= 4
        switch interaction
            case 'IV'

                d.Z2     =  d.ltariff_raw;
                [Zknots2b] = evknots(knots,d.Z2,type) ; 
                Zknots = [ Zknots1b Zknots2b];
                d.Z = Zknots;
                d.ZA_K  = [Zknots  d.FEA  ];  
                d.ZB_K  = [Zknots  d.FEB  ];  
            case 'wealth_oIV'

                d.Z2     =  d.ltariff_raw;
                [Zknots2b] = evknots(knots,d.Z2,type) ; 
                Zknots = [ Zknots1b Zknots2b];
                d.Z = Zknots;

                disp('discrete measure of origin wealth')
                d.category = d.rich_I;
                make_knots_dual_combo

                dc0 = repmat(d.category==0,1,size(Zknots,2));
                dc1 = repmat(d.category==1,1,size(Zknots,2));
                d.ZA_K =   [Zknots.*dc0  Zknots.*dc1   d.FEA ];  
                d.ZB_K =   [Zknots.*dc0  Zknots.*dc1   d.FEB ]; 

                dc0 = repmat(d.category==0,1,size(Xknots,2));
                dc1 = repmat(d.category==1,1,size(Xknots,2));
                d.X_K =   [ Xknots.*dc0  Xknots.*dc1  ];
                d.Deriv = [d.Deriv.*dc0 d.Deriv.*dc1  ];

                d.type1 = 'Developed Origin';
                d.type0 = 'Developing Origin';

            case 'comlang_ethno'
                disp('discrete measure of common language (shared by 9% or more)')
                d.category = d.comlang;
                make_knots_dual_combo
                d.ZA_K = [Zknots Zknots.*repmat(d.comlang,1,size(Zknots,2))  d.FEA  d.comlang];  
                d.ZB_K = [Zknots Zknots.*repmat(d.comlang,1,size(Zknots,2))  d.FEB  d.comlang];  
                d.X_K  = [Xknots Xknots.*repmat(d.comlang,1,size(Xknots,2)) ];
                d.Deriv = [d.Deriv  d.Deriv.*repmat(d.category,1,size(Xknots,2)) ];
                d.type1 = 'Common Language';
                d.type0 = 'No Common Language';
            case 'comcur'
                disp('discrete measure of common currency')
                d.category = d.comcur;
                make_knots_dual_combo
                d.ZA_K = [Zknots Zknots.*repmat(d.comcur,1,size(Zknots,2))  d.FEA  d.comcur];  
                d.ZB_K = [Zknots Zknots.*repmat(d.comcur,1,size(Zknots,2))  d.FEB  d.comcur];  
                d.X_K  = [Xknots Xknots.*repmat(d.comcur,1,size(Xknots,2)) ];
                d.Deriv = [d.Deriv  d.Deriv.*repmat(d.category,1,size(Xknots,2)) ];
                d.type1 = 'Common Currency';
                d.type0 = 'No Common Currency';

            case 'deep'
                disp('discrete measure of common currency + FTA')
                d.category = d.comcur & d.fta_wto;
                make_knots_dual_combo
                d.ZA_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2))  d.FEA  d.category];  
                d.ZB_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2))  d.FEB  d.category];  
                d.X_K  = [Xknots Xknots.*repmat(d.category,1,size(Xknots,2)) ];
                d.Deriv = [d.Deriv  d.Deriv.*repmat(d.category,1,size(Xknots,2)) ];
                d.type1 = 'Deep Integration';
                d.type0 = 'No Deep Integration';
            case 'langcol'
                disp('discrete measure of Language + Colonial')
                d.category = d.colony | d.comlang;
                make_knots_dual_combo
                d.ZA_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2))  d.FEA  d.category];  
                d.ZB_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2))  d.FEB  d.category];  
                d.X_K  = [Xknots Xknots.*repmat(d.category,1,size(Xknots,2)) ];
                d.Deriv = [d.Deriv  d.Deriv.*repmat(d.category,1,size(Xknots,2)) ];
                d.type1 = 'Language or Colony';
                d.type0 = 'No Language or Colony';

            case 'colony'
                disp('discrete measure of colonial relationship')
                d.category = d.colony;
                make_knots_dual_combo
                d.ZA_K = [Zknots Zknots.*repmat(d.colony,1,size(Zknots,2))  d.FEA  d.colony];  
                d.ZB_K = [Zknots Zknots.*repmat(d.colony,1,size(Zknots,2))  d.FEB  d.colony];  
                d.X_K  = [Xknots Xknots.*repmat(d.colony,1,size(Xknots,2)) ];
                d.Deriv = [d.Deriv  d.Deriv.*repmat(d.category,1,size(Xknots,2)) ];
                d.type1 = 'Colonial Relationship';
                d.type0 = 'No Colonial Relationship';
            case 'contig'
                disp('discrete measure of Border relationship')
                d.category = d.contig;
                make_knots_dual_combo
                d.FE_C = [d.contig];
                d.ZA_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2)) d.FEA  d.contig];  
                d.ZB_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2)) d.FEB  d.contig];  
                d.X_K  = [Xknots Xknots.*repmat(d.category,1,size(Xknots,2))];
                d.Deriv = [d.Deriv  d.Deriv.*repmat(d.category,1,size(Xknots,2)) ];
                d.type1 = 'Border Countries';
                d.type0 = 'Non-Border Countries';
            case 'fta_wto'
                disp('discrete measure of FTA, measured by WTO (via Penn World Tables)')
                d.category = d.fta_wto;
                make_knots_dual_combo
                d.FE_C = [d.fta_wto];
                d.ZA_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2)) d.FEA  ];  
                d.ZB_K = [Zknots Zknots.*repmat(d.category,1,size(Zknots,2)) d.FEB  ];  
                d.X_K  = [Xknots Xknots.*repmat(d.category,1,size(Xknots,2))];
                d.Deriv = [d.Deriv  d.Deriv.*repmat(d.category,1,size(Xknots,2)) ];
                d.type1 = 'FTA (defined by WTO)';
                d.type0 = 'non-FTA';
            case 'wealth_dI'
                disp('discrete measure of destination wealth')
                d.category = d.rich_J;
                make_knots_dual_combo

                dc0 = repmat(d.category==0,1,size(Zknots,2));
                dc1 = repmat(d.category==1,1,size(Zknots,2));
                d.ZA_K =   [Zknots.*dc0  Zknots.*dc1   d.FEA ];  
                d.ZB_K =   [Zknots.*dc0  Zknots.*dc1   d.FEB ]; 

                dc0 = repmat(d.category==0,1,size(Xknots,2));
                dc1 = repmat(d.category==1,1,size(Xknots,2));
                d.X_K =   [ Xknots.*dc0  Xknots.*dc1  ];
                d.Deriv = [d.Deriv.*dc0 d.Deriv.*dc1  ];

                d.type1 = 'Developed Destinations';
                d.type0 = 'Developing Destinations';
            case 'wealth_oI'
                disp('discrete measure of origin wealth')
                d.category = d.rich_I;
                make_knots_dual_combo

                dc0 = repmat(d.category==0,1,size(Zknots,2));
                dc1 = repmat(d.category==1,1,size(Zknots,2));
                d.ZA_K =   [Zknots.*dc0  Zknots.*dc1   d.FEA ];  
                d.ZB_K =   [Zknots.*dc0  Zknots.*dc1   d.FEB ]; 

                dc0 = repmat(d.category==0,1,size(Xknots,2));
                dc1 = repmat(d.category==1,1,size(Xknots,2));
                d.X_K =   [ Xknots.*dc0  Xknots.*dc1  ];
                d.Deriv = [d.Deriv.*dc0 d.Deriv.*dc1  ];

                d.type1 = 'Developed Origin';
                d.type0 = 'Developing Origin';
            case 'china_oI'
                disp('discrete measure of China (1/0)')
                d.category = d.rich_I;
                make_knots_dual_combo

                dc0 = repmat(d.category==0,1,size(Zknots,2));
                dc1 = repmat(d.category==1,1,size(Zknots,2));
                d.ZA_K =   [Zknots.*dc0  Zknots.*dc1   d.FEA ];  
                d.ZB_K =   [Zknots.*dc0  Zknots.*dc1   d.FEB ]; 

                dc0 = repmat(d.category==0,1,size(Xknots,2));
                dc1 = repmat(d.category==1,1,size(Xknots,2));
                d.X_K =   [ Xknots.*dc0  Xknots.*dc1  ];
                d.Deriv = [d.Deriv.*dc0 d.Deriv.*dc1  ];

                d.type1 = 'China Origin';
                d.type0 = 'Other Origin';
            case 'wealth_odI'
                disp('discrete measure of orign and destination wealth')
                d.category = zeros(size(d.rich_I));
                d.category(d.rich_I == 1 & d.rich_J == 1) = 0;
                d.category(d.rich_I == 1 & d.rich_J == 0) = 1;
                d.category(d.rich_I == 0 & d.rich_J == 1) = 2;
                d.category(d.rich_I == 0 & d.rich_J == 0) = 3;

                d.type0 = 'Developed -> Developed';
                d.type1 = 'Developed -> Developing';
                d.type2 = 'Developing -> Developed';
                d.type3 = 'Developing -> Developing';

                make_knots_quad_combo

                dc0 = repmat(d.category==0,1,size(Zknots,2));
                dc1 = repmat(d.category==1,1,size(Zknots,2));
                dc2 = repmat(d.category==2,1,size(Zknots,2));
                dc3 = repmat(d.category==3,1,size(Zknots,2));
              
                d.ZA_K =   [Zknots.*dc0  Zknots.*dc1  Zknots.*dc2  Zknots.*dc3  d.FEA  d.FE_C];  
                d.ZB_K =   [Zknots.*dc0  Zknots.*dc1  Zknots.*dc2  Zknots.*dc3  d.FEB  d.FE_C];  

                dc0 = repmat(d.category==0,1,size(Xknots,2));
                dc1 = repmat(d.category==1,1,size(Xknots,2));
                dc2 = repmat(d.category==2,1,size(Xknots,2));
                dc3 = repmat(d.category==3,1,size(Xknots,2));

                d.X_K =   [ Xknots.*dc0  Xknots.*dc1  Xknots.*dc2  Xknots.*dc3];
                d.Deriv = [d.Deriv.*dc0 d.Deriv.*dc1 d.Deriv.*dc2 d.Deriv.*dc3];
            case 'wealth_ooI'
                disp('discrete measure of orign wealth')



                % World Bank Cutoffs for 2000
                %                 <= 755
                %                 756-2,995
                %                 2,996-9,265
                %                 > 9,265
                d.cutoff1 = log(755);
                d.cutoff2 = log(2995);
                d.cutoff3 = log(9265);


                d.category = zeros(size(d.rich_I));
                d.category(d.rich_I == 1 ) = 3;
                d.category(d.rich_I == 0 & d.lGDPc_I >  d.cutoff2) = 2;
                d.category(d.rich_I == 0 & d.lGDPc_I >  d.cutoff1 & d.lGDPc_I <= d.cutoff2) = 1;
                d.category(d.rich_I == 0 & d.lGDPc_I <= d.cutoff1) = 0;
                hist(d.category)

                d.type3 = 'High';
                d.type2 = 'Upper middle';
                d.type1 = 'Lower middle';
                d.type0 = 'Low';

                make_knots_quad_combo

                dc0 = repmat(d.category==0,1,size(Zknots,2));
                dc1 = repmat(d.category==1,1,size(Zknots,2));
                dc2 = repmat(d.category==2,1,size(Zknots,2));
                dc3 = repmat(d.category==3,1,size(Zknots,2));
              
                d.ZA_K =   [Zknots.*dc0  Zknots.*dc1  Zknots.*dc2  Zknots.*dc3  d.FEA  d.FE_C];  
                d.ZB_K =   [Zknots.*dc0  Zknots.*dc1  Zknots.*dc2  Zknots.*dc3  d.FEB  d.FE_C];  

                dc0 = repmat(d.category==0,1,size(Xknots,2));
                dc1 = repmat(d.category==1,1,size(Xknots,2));
                dc2 = repmat(d.category==2,1,size(Xknots,2));
                dc3 = repmat(d.category==3,1,size(Xknots,2));

                d.X_K =   [ Xknots.*dc0  Xknots.*dc1  Xknots.*dc2  Xknots.*dc3];
                d.Deriv = [d.Deriv.*dc0 d.Deriv.*dc1 d.Deriv.*dc2 d.Deriv.*dc3];

            otherwise
                disp('non-valid interaction')
        end
end
    


%% Pre-Compute Projection Matrices (to save on computing cost down the road)
    disp('create projection')
    d.Z_K       = [d.ZA_K;d.ZB_K];
    d.PZ_K      = d.Z_K*(d.Z_K'*d.Z_K)^-1*d.Z_K';
    d.FE        = [[d.FEA d.FE_C];[d.FEB d.FE_C]];
    d.FE_PZ_Ki  = (d.FE'*d.PZ_K*d.FE)^-1*d.FE'*d.PZ_K;

%     CHECK IF THIS SEEMS REASONABLE (SHOULDN'T BE TOO BIG)
    sum(sum(inv(d.Z_K'*d.Z_K)))


end


