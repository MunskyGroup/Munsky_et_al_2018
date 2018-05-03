function Update_Plots
% Compute and Save Means, Moments, Distribution Analyses, Spatial and
% NonSpatial, for STL1 and CTT1.
%%
for FT = {'Moments_4','Moments_Corrected','Distributions','Means'}
% for FT = {'Moments_Corrected','Distributions','Means'}
    switch FT{1}
        case {'Moments_4','Moments_Corrected','Means'}
%             GENES_AR = {'STL1_and_CTT1','CTT1','STL1'};
            GENES_AR = {'CTT1','STL1'};
        case 'Distributions'
            GENES_AR = {'CTT1','STL1'};
    end
    for SP={'NonSpatial','Spatial'}
        for GENES = GENES_AR
            close all
            if strcmp(GENES{1},'STL1_and_CTT1')
                Model = '4State_2gene_dep_k21only_noloop';
                DIR = ['Fit_Results/',GENES{1},'/',SP{1},'/',FT{1},'/',Model];
            else
                Model = '4State_k21only_noloop';
                DIR = ['Fit_Results/',GENES{1},'_0p2and0p4/',SP{1},'/',FT{1},'/',Model];
            end
            
            Results_File = [DIR,'/Analyses.mat'];
            load(Results_File)
            Hog.MOD.Type='Moments_Corrected';
            Hog_Model.plot_moments_mod(Hog.MOD,Hog.Moments,1,2);
            for TMP = {'_A','_B','_C'}
                Hog.MOD.Reps = TMP{1};
                Hog_Model.plot_moments_dat(Hog.MOD,Hog.MOD.Data,1,2);
            end
            
            f1=figure(1); f2=figure(2);
            if strcmp(GENES{1},'STL1_and_CTT1')
                f1 = Adjust_Plots.Mean_Joint(f1,Hog.MOD);
                f2 = Adjust_Plots.Vars_Joint(f2,Hog.MOD);
            else
                Adjust_Plots.Mean(1,Hog.MOD);
                Adjust_Plots.Vars(2,Hog.MOD);
            end
            
            if ~strcmp(GENES{1},'STL1_and_CTT1')  % Distributions found with FSP only for the single gene models.
                Hog.MOD.Type='Distributions';
                Hog_Model.plot_distributions_mod(Hog.MOD,Hog.Distributions,3,5);
                for TMP = {'_A','_B','_C'}
                    Hog.MOD.Reps = TMP{1};
                    Hog_Model.plot_distributions_dat(Hog.MOD,Hog.MOD.Data,3,5);
                end
                Adjust_Plots.Marginal(3,Hog.MOD,0.2)
                Adjust_Plots.Marginal(4,Hog.MOD,0.4)
                Adjust_Plots.On(5,'NonSpatial')
                
                if strcmp(SP{1},'Spatial')
                    Hog.MOD.Type='Distributions';
                    Hog_Model.plot_spat_distributions_mod(Hog.MOD,Hog.Distributions,7,11);
                    for TMP = {'_A','_B','_C'}
                        Hog.MOD.Reps = TMP{1};
                        Hog_Model.plot_spat_distributions_dat(Hog.MOD,Hog.MOD.Data,9,11);
                    end
                    Adjust_Plots.On(11,'Spatial')
                    
                    Hog.MOD.Reps = '_AC';
                    Hog_Model.plot_spat_distributions_dat(Hog.MOD,Hog.MOD.Data,9,100); close(100);
                    Adjust_Plots.Joint(7,Hog.MOD,0.2,'Model','Cyt','Nuc')
                    Adjust_Plots.Joint(8,Hog.MOD,0.4,'Model','Cyt','Nuc')
                    Adjust_Plots.Joint(9,Hog.MOD,0.2,'Data','Cyt','Nuc')
                    Adjust_Plots.Joint(10,Hog.MOD,0.4,'Data','Cyt','Nuc')
                end
            elseif  strcmp(GENES{1},'STL1_and_CTT1') % Use SSA results.
                switch SP{1}
                    case 'NonSpatial'
                        Hog_Model.plot_SSA_dist_results(Hog.MOD,Hog.SSA,13,5,[1],[])
                        obj = Hog.MOD;
                        obj.Type='Distributions';
                        obj.Gene = 'STL1';
                        for TMP = {'_A','_B','_C'}
                            obj.Reps = TMP{1};
                            Hog_Model.plot_distributions_dat(obj,obj.Data,13,5);
                        end
                        Adjust_Plots.Marginal(13,Hog.MOD,0.2)
                        Adjust_Plots.Marginal(14,Hog.MOD,0.2)
                        Adjust_Plots.On(5,'NonSpatial')
                        
                        Hog_Model.plot_SSA_dist_results(Hog.MOD,Hog.SSA,23,5,[2],[])
                        obj.Gene = 'CTT1';
                        for TMP = {'_A','_B','_C'}
                            obj.Reps = TMP{1};
                            Hog_Model.plot_distributions_dat(obj,obj.Data,23,5);
                        end
                        Adjust_Plots.Marginal(23,Hog.MOD,0.2)
                        Adjust_Plots.Marginal(24,Hog.MOD,0.2)
                        Adjust_Plots.On(5,'NonSpatial')

                    case {'Spatial','Fixed'}
                        Hog_Model.plot_SSA_dist_results(Hog.MOD,Hog.SSA,13,5,[1,3],[])
                        obj = Hog.MOD;
                        obj.Type='Distributions';
                        obj.Gene = 'STL1';
                        for TMP = {'_A','_B','_C'}
                            obj.Reps = TMP{1};
                            Hog_Model.plot_distributions_dat(obj,obj.Data,13,5);
                        end
                        Adjust_Plots.Marginal(13,Hog.MOD,0.2)
                        Adjust_Plots.Marginal(14,Hog.MOD,0.2)
                        Adjust_Plots.On(5,'NonSpatial')
                        
                        Hog_Model.plot_SSA_dist_results(Hog.MOD,Hog.SSA,23,5,[2,4],[])
                        obj.Gene = 'CTT1';
                        for TMP = {'_A','_B','_C'}
                            obj.Reps = TMP{1};
                            Hog_Model.plot_distributions_dat(obj,obj.Data,23,5);
                        end
                        Adjust_Plots.Marginal(23,Hog.MOD,0.2)
                        Adjust_Plots.Marginal(24,Hog.MOD,0.2)
                        Adjust_Plots.On(5,'NonSpatial')

                end
            end
            
            FIG_Sub = ['Figures/',GENES{1},'_',SP{1},'_Fit_',FT{1},'_Plot_'];
            
            saveas(f1,[FIG_Sub,'Means.fig']);
            saveas(f2,[FIG_Sub,'Vars.fig']);
            try saveas(5,[FIG_Sub,'On_Total.fig']); catch; end;
            try saveas(11,[FIG_Sub,'On_Spatial.fig']); catch; end;
            try saveas(3,[FIG_Sub,'Marg_Dists_0p2.fig']); catch; end;
            try saveas(4,[FIG_Sub,'Marg_Dists_0p4.fig']); catch; end;
            try saveas(13,[FIG_Sub,'Marg_Dists_0p2_STL1.fig']); catch; end;
            try saveas(14,[FIG_Sub,'Marg_Dists_0p4_STL1.fig']); catch; end;
            try saveas(23,[FIG_Sub,'Marg_Dists_0p2_CTT1.fig']); catch; end;
            try saveas(24,[FIG_Sub,'Marg_Dists_0p4_CTT1.fig']); catch; end;
            try saveas(7,[FIG_Sub,'Spat_Dists_Mod_0p2.fig']); catch; end;
            try saveas(8,[FIG_Sub,'Spat_Dists_Mod_0p4.fig']); catch; end;
            try saveas(9,['Figures/',GENES{1},'_Spat_Dists_Dat_0p2.fig']); catch; end;
            try saveas(10,['Figures/',GENES{1},'_Spat_Dists_Dat_0p4.fig']); catch; end;
        end
    end
end

