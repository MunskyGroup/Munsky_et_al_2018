
close all
if ~exist('RNAdata','var')
    load('/Users/munsky/Dropbox/Shared_Codes/Hog_Model/Hog_Data/RNAspots.mat')
end

for type_ar = {'Distributions'}
    
    type = type_ar{1};
    fil = ['Fit_Results/K_elong_rates_',type,'.mat'];
    %% First for all data and both genes.
%     Elong_All = get_TS_Fit_errs_All(RNAdata,'BOTH',[1 2],[1 2 3],'Spatial',type);
%     save(fil,'Elo*');
    
    %% First for all data
    Elong_CTT1 = get_TS_Fit_errs_All(RNAdata,'CTT1',[1 2],[1 2 3],'Spatial',type);
    Elong_STL1 = get_TS_Fit_errs_All(RNAdata,'STL1',[1 2],[1 2 3],'Spatial',type);
    save(fil,'Elo*');
    
    %% Next for all data at each salt level.
    Elong_CTT1_02 = get_TS_Fit_errs_All(RNAdata,'CTT1',[1],[1 2],'Spatial',type);
    Elong_STL1_02 = get_TS_Fit_errs_All(RNAdata,'STL1',[1],[1 2],'Spatial',type);
    
    Elong_CTT1_04 = get_TS_Fit_errs_All(RNAdata,'CTT1',[2],[1 2 3],'Spatial',type);
    Elong_STL1_04 = get_TS_Fit_errs_All(RNAdata,'STL1',[2],[1 2 3],'Spatial',type);
    save(fil,'Elo*');
    
    %% Next, for every individual combination.
    clear Elongs
    for l=1:2
        switch l
            case 1
                spat='Spatial';
            case 2
                spat='NonSpatial';
        end
        for k=1:2
            switch k
                case 1
                    g='STL1';
                case 2
                    g='CTT1';
            end
            for i=1:2
                switch i
                    case 1
                        nreps = 2;
                    case 2
                        nreps = 3;
                end
                for j=1:nreps
                    a = get_TS_Fit_errs_All(RNAdata,g,i,j,spat,type);
                    Elongs(k,i,j,l) = a;
                end
            end
        end
    end
    
    save(fil,'Elo*');
    
end