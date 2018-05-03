function MCMC_Chain_Obj = Process_MCMC_class(obj)

%% Load the MCMC file if it is precompiled.
if (obj.Recompile_MCMC==0)&&~isempty(obj.MCMC_File)&&exist(obj.MCMC_File,'file')
    load(obj.MCMC_File,'MCMC_Chain_Obj');
    
elseif obj.Recompile_MCMC==0
    error('Cannot find the appropriate processed MCMC File')
    
elseif obj.Recompile_MCMC==1
    switch obj.Model
        case '4State_k21only_noloop'
            switch obj.Spatial
                case 'Spatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k21b','kr1','kr2','kr3','kr4','\gamma_n','t0','\gamma_c','tpt'};
                    MCMC_Chain_Obj.Pos_Pars = [1,3,4,5,6,8:15];
                case 'Fixed'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k21b','kr1','kr2','kr3','kr4','\gamma_n','t0','tpt'};
                    MCMC_Chain_Obj.Pos_Pars = [1,3,4,5,6,8:14];
                case 'NonSpatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k21b','kr1','kr2','kr3','kr4','\gamma','t0'};
                    MCMC_Chain_Obj.Pos_Pars = [1,3,4,5,6,8:13];
            end
        case '4State_Full'
            switch obj.Spatial
                case 'Spatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k12b','k21b','k23b','k32b','k34b','k43b','kr1','kr2','kr3','kr4','\gamma_n','t0','\gamma_c','tpt',};
                    MCMC_Chain_Obj.Pos_Pars = [13:20];
                case 'Fixed'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k12b','k21b','k23b','k32b','k34b','k43b','kr1','kr2','kr3','kr4','\gamma','t0','tpt',};
                    MCMC_Chain_Obj.Pos_Pars = [13:19];
                case 'NonSpatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k12b','k21b','k23b','k32b','k34b','k43b','kr1','kr2','kr3','kr4','\gamma','t0'};
                    MCMC_Chain_Obj.Pos_Pars = [13:18];
            end
        case '3State_Full'
            switch obj.Spatial
                case 'Spatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k13a','k21a','k23a','k31a','k32a','k12b','k13b','k21b','k23b','k31b','k32b','kr1','kr2','kr3','\gamma_n','t0','\gamma_c','tpt',};
                    MCMC_Chain_Obj.Pos_Pars = [13:19];
                case 'Fixed'
                    MCMC_Chain_Obj.par_names = {'k12a','k13a','k21a','k23a','k31a','k32a','k12b','k13b','k21b','k23b','k31b','k32b','kr1','kr2','kr3','\gamma','t0','tpt',};
                    MCMC_Chain_Obj.Pos_Pars = [13:18];
                case 'NonSpatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k13a','k21a','k23a','k31a','k32a','k12b','k13b','k21b','k23b','k31b','k32b','kr1','kr2','kr3','\gamma','t0'};
                    MCMC_Chain_Obj.Pos_Pars = [13:17];
            end
        case '4State_2gene_dep_k21only_noloop'
            switch obj.Spatial
                case 'Spatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k21b','kr1','kr2','kr3','kr4','\gamma_n','t0','\gamma_c','tpt'};
                    MCMC_Chain_Obj.Pos_Pars = [1,3,4,5,6,8:15];
                case 'Fixed'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k21b','kr1','kr2','kr3','kr4','\gamma_n','t0','tpt'};
                    MCMC_Chain_Obj.Pos_Pars = [1,3,4,5,6,8:14];
                case 'NonSpatial'
                    MCMC_Chain_Obj.par_names = {'k12a','k21a','k23a','k32a','k34a','k43a','k21b','kr1','kr2','kr3','kr4','\gamma','t0'};
                    MCMC_Chain_Obj.Pos_Pars = [1,3,4,5,6,8:13];
            end
    end
    
    switch obj.Fit_Type
        case 'Means'
            switch obj.Gene
                case 'STL1'
                    switch obj.Spatial
                        case 'Spatial'
                            INDS = [211:216,218,219,220];
                        case 'NonSpatial'
                            INDS = [202,204:210];
                    end
                case 'CTT1'
                    switch obj.Spatial
                        case 'Spatial'
                            INDS = [221:230];
                        case 'NonSpatial'
                            INDS = [202:210];
                    end
            end
        case 'Moments_Corrected'
            switch obj.Gene
                case 'STL1'
                    switch obj.Spatial
                        case 'Spatial'
                            INDS = [211:220];
                        case 'NonSpatial'
                            INDS = [222:230];
                    end
                case 'CTT1'
                    switch obj.Spatial
                        case 'Spatial'
                            INDS = [221:230];
                        case 'NonSpatial'
                            INDS = [202:210];
                    end
            end
        case 'Moments_4'
            switch obj.Spatial
                case 'Spatial'
                    switch obj.Gene
                        case 'STL1'
                            INDS = [211:220];
                        case 'CTT1'
                            INDS = [221:230];
                    end
                case 'NonSpatial'
                    switch obj.Gene
                        case 'STL1'
                            INDS = [203,204,206,208,209];
                        case 'CTT1'
                            INDS = [231:240];
                    end
            end
        case 'Distributions'
            switch obj.Spatial
                case 'Spatial'
                    switch obj.Gene
                        case 'STL1'
                            INDS = [202:220];
                        case 'CTT1'
                            INDS = [202:220];
                    end
                case 'NonSpatial'
                    switch obj.Gene
                        case 'STL1'
                            INDS = [202:207,209,210];
                        case 'CTT1'
                            INDS = [202:220];
                    end
            end
        case 'Gaussian'
            INDS = [201:210];
    end
    skip_pt = 10;
    BestVal=-inf;
    for i_ind=1:length(INDS)
        i = INDS(i_ind);
        biggest = [];
        try
            fnms = ls([obj.MCMC_Raw_Dir,'/Hog_Met_Hast_N_',num2str(i),'_pt_*']);
            
            for j=6000:-1:1
                File_Name = [obj.MCMC_Raw_Dir,'/Hog_Met_Hast_N_',num2str(i),'_pt_',num2str(j),obj.Reps,'.mat'];
                if contains(fnms,File_Name)
                    biggest = j;
                    break
                end
            end
        catch
        end
        if ~isempty(biggest)
            Met_Hast_Sample=[];
            File_Name = [obj.MCMC_Raw_Dir,'/Hog_Met_Hast_N_',num2str(i),'_pt_',num2str(biggest),obj.Reps,'.mat'];
            load(File_Name)
            bn = biggest*500;
            
            Samples = zeros(length(bn:skip_pt:biggest*1000),length(SCALE));
            for i1=1:length(SCALE)
                DATX = Met_Hast_Sample(bn:skip_pt:biggest*1000,i1)*SCALE(i1);
                if min((i1-MCMC_Chain_Obj.Pos_Pars).^2)==0
                    DATX = 10.^DATX;
                end
                Samples(:,i1) = DATX;
            end
            
            Met_Hast_Sample_Values = Met_Hast_Sample_Values(1:biggest*1000,1);
            MCMC_Chain_Obj.chain(i_ind).obj = Met_Hast_Sample_Values(bn:skip_pt:biggest*1000,1);
            MCMC_Chain_Obj.chain(i_ind).pars = Samples;
            
            [BB,J] = max(Met_Hast_Sample_Values(isfinite(Met_Hast_Sample_Values)));
            
            if BB>BestVal
                BestVal=BB;
                Best = Met_Hast_Sample(J,:);
                for i1=1:length(SCALE)
                    Best(i1) = Best(i1)*SCALE(i1);
                    if min((i1-MCMC_Chain_Obj.Pos_Pars).^2)==0
                        Best(i1) = 10.^Best(i1);
                    end
                end
                MCMC_Chain_Obj.BestVal = BestVal;
                MCMC_Chain_Obj.Best = Best;
            end
        end
    end
    save(obj.MCMC_File,'MCMC_Chain_Obj');
end



