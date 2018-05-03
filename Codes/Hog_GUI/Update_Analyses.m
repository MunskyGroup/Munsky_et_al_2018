function Update_Analyses
% Compute and Save Means, Moments, Distribution Analyses, Spatial and
% NonSpatial, for STL1 and CTT1.

MOD = Hog_Model;
MOD.Salt = [0.2 0.4];
parfor i=1:20
    Save_Rslts(MOD,i)
end
Make_Table

function Save_Rslts(MOD,i)
I(1) = floor((i-1)/4)+1;
b = mod(i-1,4)+1;
I(2) = floor((b-1)/2)+1;
I(3) = mod(b-1,2)+1;
Hog.MOD = MOD;

FT = {'Gaussian','Moments_4','Moments_Corrected','Means','Distributions'};
Hog.Fit_Type = FT{I(1)};

SP={'NonSpatial','Spatial'};
Hog.MOD.Spatial=SP{I(2)};

GENES = {'CTT1','STL1'};%,'STL1_and_CTT1'};
Hog.MOD.Gene = GENES{I(3)};

switch GENES{I(3)}
    case {'CTT1','STL1'}
        Hog.MOD.Model = '4State_k21only_noloop';
        Hog.MOD.MGenes = 1;
        DIR = ['Fit_Results/',Hog.MOD.Gene,'_0p2and0p4/',Hog.MOD.Spatial,'/',Hog.Fit_Type,'/',Hog.MOD.Model];
    case 'STL1_and_CTT1'
        Hog.MOD.Model = '4State_2gene_dep_k21only_noloop';
        Hog.MOD.MGenes = 2;
        DIR = ['Fit_Results/',Hog.MOD.Gene,'/',Hog.MOD.Spatial,'/',Hog.Fit_Type,'/',Hog.MOD.Model];
end

Hog.MOD.Parameter_File = [DIR,'/Best_Pars_AC.mat'];

if strcmp(Hog.Fit_Type,'Distributions')
    disp('Starting some SSA analyses needed for Figs. 3 and S15. This could take a while.');
    Hog.MOD.Type = 'SSA';
    Hog.MOD.SSA_Runs=10000;
    Hog.MOD.SSA_Cells = 1000;
    Hog.SSA = Hog_Model.solve_hog(Hog.MOD);
end

Hog.MOD.Type='Moments';
Hog.MOD.Fit_Type='Moments';
Hog.Moments = Hog_Model.solve_hog(Hog.MOD);

Hog.MOD.Fit_Type='Moments_4';
Hog.MOD.Type='Moments_4';
Hog.MOD.Moment_Order = 4;
Hog.Moments_4 = Hog_Model.solve_hog(Hog.MOD);

if ~strcmp(GENES{I(3)},'STL1_and_CTT1')
    Hog.MOD.Fit_Type='Distributions';
    Hog.MOD.Type='Distributions';
    Hog.Distributions = Hog_Model.solve_hog(Hog.MOD);
end

Results_File = [DIR,'/Analyses.mat'];
save(Results_File,'Hog','-v7.3');
disp(['Completed ', DIR]);

function Make_Table
FT = {'Means','Moments_Corrected','Moments_4','Distributions'};

SP={'NonSpatial','Spatial'};

GENES = {'CTT1','STL1'};%,'STL1_and_CTT1'};


for j=1:2
    for k=1:2
        for i = 1:4
            Hog.Fit_Type = FT{i};
            Hog.MOD.Spatial=SP{k};
            Hog.MOD.Gene = GENES{j};
            
            switch Hog.MOD.Gene
                case {'CTT1','STL1'}
                    Hog.MOD.Model = '4State_k21only_noloop';
                    Hog.MOD.MGenes = 1;
                    DIR = ['Fit_Results/',Hog.MOD.Gene,'_0p2and0p4/',Hog.MOD.Spatial,'/',Hog.Fit_Type,'/',Hog.MOD.Model];
                case 'STL1_and_CTT1'
                    Hog.MOD.Model = '4State_2gene_dep_k21only_noloop';
                    Hog.MOD.MGenes = 2;
                    DIR = ['Fit_Results/',Hog.MOD.Gene,'/',Hog.MOD.Spatial,'/',Hog.Fit_Type,'/',Hog.MOD.Model];
            end
            Results_File = [DIR,'/Analyses.mat']
            
            load(Results_File,'Hog');
            
            
            AAA(j,k,i,1) = Hog.Moments(1).likelihoodC+Hog.Moments(2).likelihoodC;
            AAA(j,k,i,2) = Hog.Moments(1).likelihood+Hog.Moments(2).likelihood;
            AAA(j,k,i,3) = Hog.Moments_4(1).likelihood+Hog.Moments_4(2).likelihood;
            AAA(j,k,i,4) = Hog.Distributions(1).likelihood+Hog.Distributions(2).likelihood;
            
        end
    end
end

fid = fopen(['Likelihoods.csv'],'w');
fprintf(fid,'gene,&,analysis,&,mn,&,mn and var,&,Lygeros,&,dist');
for j=1:2
    for k=1:2
        for i = 1:4
            Hog.Fit_Type = FT{i};
            Hog.MOD.Spatial=SP{k};
            Hog.MOD.Gene = GENES{j};
            
            switch Hog.MOD.Gene
                case {'CTT1','STL1'}
                    Hog.MOD.Model = '4State_k21only_noloop';
                    Hog.MOD.MGenes = 1;
                    DIR = ['Fit_Results/',Hog.MOD.Gene,'_0p2and0p4/',Hog.MOD.Spatial,'/',Hog.Fit_Type,'/',Hog.MOD.Model];
                case 'STL1_and_CTT1'
                    Hog.MOD.Model = '4State_2gene_dep_k21only_noloop';
                    Hog.MOD.MGenes = 2;
                    DIR = ['Fit_Results/',Hog.MOD.Gene,'/',Hog.MOD.Spatial,'/',Hog.Fit_Type,'/',Hog.MOD.Model];
            end
            Results_File = [DIR,'/Analyses.mat']
            
            load(Results_File,'Hog');
            fprintf(fid,['\n ',GENES{j},',&,',FT{i},',&,',SP{k},',&,']);
            fprintf(fid,num2str(-Hog.Moments(1).likelihoodC-Hog.Moments(2).likelihoodC+min(squeeze(AAA(j,k,:,1)))));fprintf(fid,',&,');
            fprintf(fid,num2str(-Hog.Moments(1).likelihood-Hog.Moments(2).likelihood+min(squeeze(AAA(j,k,:,2)))));fprintf(fid,',&,');
            fprintf(fid,num2str(-Hog.Moments_4(1).likelihood-Hog.Moments_4(2).likelihood+min(squeeze(AAA(j,k,:,3)))));fprintf(fid,',&,');
            %             fprintf(fid,num2str(Hog.Moments.likelihoodC));fprintf(fid,',&,');
            fprintf(fid,num2str(-Hog.Distributions(1).likelihood-Hog.Distributions(2).likelihood+min(squeeze(AAA(j,k,:,4)))));
            
        end
    end
end

