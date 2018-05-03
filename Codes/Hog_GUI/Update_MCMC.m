MOD = Hog_Model;
MOD.Recompile_MCMC = 1;
for FT = {'Gaussian','Distributions','Means','Moments_Corrected','Moments_4',}
% for FT = {'Distributions'}
% for FT = {'Distributions','Means','Moments_Corrected','Moments_4'}
    Hog.MOD = MOD;
    Hog.MOD.Fit_Type = FT{1};
    
    for SP={'NonSpatial','Spatial'}
        Hog.MOD.Spatial=SP{1};
        
        for GENES = {'CTT1','STL1'};%,'STL1_and_CTT1'};
            Hog.MOD.Gene = GENES{1};
            if strcmp(GENES{1},'STL1_and_CTT1')
                Hog.MOD.Model = '4State_2gene_dep_k21only_noloop';
                Hog.MOD.MCMC_Raw_Dir = ['/Users/munsky/Dropbox/Shared_Codes/SSIT_Data_and_Results/Hog_System/MCMC_Results/',GENES{1},'/',SP{1},'/',FT{1},'/',Hog.MOD.Model];                
                DIR = ['Fit_Results/',Hog.MOD.Gene,'/',Hog.MOD.Spatial,'/',Hog.MOD.Fit_Type,'/',Hog.MOD.Model];
            else
                Hog.MOD.Model ='4State_k21only_noloop';
                Hog.MOD.MCMC_Raw_Dir = ['/Users/munsky/Dropbox/Shared_Codes/SSIT_Data_and_Results/Hog_System/MCMC_Results/',GENES{1},'_0p2and0p4/',SP{1},'/',FT{1},'/',Hog.MOD.Model];
                Model = '4State_k21only_noloop';
                DIR = ['Fit_Results/',Hog.MOD.Gene,'_0p2and0p4/',Hog.MOD.Spatial,'/',Hog.MOD.Fit_Type,'/',Hog.MOD.Model];
            end
            Hog.MOD.MCMC_File = [DIR,'/MCMC_Processed.mat'];
            Hog.MOD.MCMC_Chain_Obj;
            
            disp(['Completed ', DIR]);

        end
    end
end