% To find the "standard" mRNa intensity, I look at all cells at all time
% points. I will then rank them from the dimmest to the brightest.  Then I
% will remove the 'n' brightest spots (in case these are TS).  Finally, I
% will find the median value of the remaining spots.
% clear all
close all
if ~exist('RNAdata','var')
    load('/Users/munsky/Dropbox/Shared_Codes/Hog_Model/Hog_Data/RNAspots.mat')
end

load '/Users/munsky/Dropbox/Shared_Codes/Hog_Model/Hog_Data/Analyzed_FISH_Data_2016_06';

load('Hog_Data/Hog_Probes.mat')

n=3;
nmax = inf;
% Number of mRNA spots required to positively ID as an mRNA. In this
% analysis, I am only going to look at cells that have more than 'n' spots
% and then use these to determine what the average mature mRNA looks like.
% This is meant to avoid using the TS when trying to determine what is a
% mature mRNA.
N = 2; % Minimum number of nascent mRNA required to be considered as a TS.
Sp_Type = 'Spatial';
speed = 63; % CTT1  - mn of on TS.
clear TMP TS_All
for ig=1:2
    switch ig
        case 1
            gene = 'STL1';
            PD = PD_STL1;
            len = 1630;
        case 2
            gene = 'CTT1';
            PD = PD_CTT1;
            len  = 1618;
    end
    
    for i_nacl = 1:2
        TS_All(i_nacl).Data = [];
        fil = ['Fit_Results/',gene,'_0p2and0p4/',Sp_Type,'/Distributions/4State_k21only_noloop/Analyses.mat'];
        try
            load(fil)
        catch
            disp('WARNING - Missing data file named:')
            fil
            if strcmp(fil,'Fit_Results/CTT1_0p2and0p4/Spatial/Distributions/4State_k21only_noloop/Analyses.mat')
                disp('You can download this file at:');
                disp('https://www.dropbox.com/s/g4sygvsh8piq0jp/Analyses.mat?dl=0')
                disp('and move it to the directory:')
                disp('Fit_Results/CTT1_0p2and0p4/Spatial/Distributions/4State_k21only_noloop/')
                disp('or you can regenerate this and all file susing the button at the top-right.')
                msgbox({'WARNING - Missing data file named',...
                    'You can download this file at:',...
                    'https://www.dropbox.com/s/g4sygvsh8piq0jp/Analyses.mat?dl=0',...
                    'and move it to the directory:',...
                    'Fit_Results/CTT1_0p2and0p4/Spatial/Distributions/4State_k21only_noloop/',...
                    'or you can regenerate this and all file susing the button at the top-right.'})
                
            else
                disp('You can regenerate this and all file susing the button at the top-right')
            end
            return
       end
        MOD_FSP = Hog.MOD;
        switch i_nacl
            case 1
                n_reps = 2;
                nacl = '0p2NaCl';
                MOD_FSP.Salt = 0.2;
                dat_col = [0.8 0 0.8];
                mod_col = [0 0 0];
                fil = ['Fit_Results/SSATS_',gene,'_Distributions_0.2_',Sp_Type,'.mat'];
            case 2
                n_reps = 3;
                nacl = '0p4NaCl';
                dat_col = [0 0.8 0];
                MOD_FSP.Salt = 0.4;
                mod_col = [0 0 0.4];
                fil = ['Fit_Results/SSATS_',gene,'_Distributions_0.4_',Sp_Type,'.mat'];
        end
        
        %% Compute FSP analysis of TS intensity
        MOD_FSP.Fit_Type='Distributions_TS';
        MOD_FSP.Type = 'Distributions_TS';
        MOD_FSP.k_elong = speed;
        MOD_FSP.gene_length=len;
        MOD_FSP.TS_Lim = 100;
        Distributions_TS = Hog_Model.solve_hog(MOD_FSP);
        load('Part_Intens.mat','Partial_Intens_Map');
        C = ones(1,N*10+1);
        Ci = ones(1,10);
        for j=N+1:MOD_FSP.TS_Lim
            C = sparse(blkdiag(C,Ci));
        end
        Partial_Intens_Map_Red = C*Partial_Intens_Map(1:MOD_FSP.TS_Lim*10+1,1:MOD_FSP.TS_Lim+1);
        inten_dist = Partial_Intens_Map_Red*Distributions_TS.distributions(1:MOD_FSP.TS_Lim+1,:);
        
        %%
        TS_all_reps=[];
        f2=figure(100*(ig-1)+10*(i_nacl-1)+2);%subplot(3,1,1);
        f3=figure(100*(ig-1)+10*(i_nacl-1)+3);
        f4=figure(100*(ig-1)+14+i_nacl);
        N_cells = 1000*ones(1,16);
        Means_Active_TS_Mod(ig,i_nacl)=Model_TS_Results([],MOD_FSP,[],PD,N,f2,f3,f4,mod_col,max(N_cells),N_cells,inten_dist,1212);
        
        load(fil,'MOD','SSA')
        for i=1:5
            Model_TS_Results_SSA(SSA,MOD,speed,PD,N,f2,f3,f4,mod_col,max(N_cells),N_cells,inten_dist);
        end
    end
end

%%
for i=[15 16 115 116]
    set(i,'position',[1245         809         274         161]);
    figure(i);set(gca,'fontsize',14,'yscale','log','ytick',10.^[-4:0],'ylim',[1e-3,1],'xlim',[0,30])
    xlabel('')
    ylabel('')
end

for i=[2 12 102 112]
    set(i,'position',[1245         809         274         161]);
    figure(i);set(gca,'fontsize',14,'yscale','linear','xlim',[0,60])
    ylabel('')
    xlabel('')
end




