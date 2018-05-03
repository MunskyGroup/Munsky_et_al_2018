% To find the "standard" mRNa intensity, I look at all cells at all time
% points. I will then rank them from the dimmest to the brightest.  Then I
% will remove the 'n' brightest spots (in case these are TS).  Finally, I
% will find the median value of the remaining spots.
if ~exist('RNAdata','var')
    load('Hog_Data/RNAspots.mat')
end

load 'Hog_Data/Analyzed_FISH_Data_2016_06';

load('Hog_Data/Hog_Probes.mat')

n=3;
nmax = inf;
% Number of mRNA spots required to positively ID as an mRNA. In this
% analysis, I am only going to look at cells that have more than 'n' spots
% and then use these to determine what the average mature mRNA looks like.
% This is meant to avoid using the TS when trying to determine what is a
% mature mRNA.
N = 2; % Minimum number of nascent mRNA required to be considered as a TS.
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
        close all
        TS_All(i_nacl).Data = [];
        fil = ['Fit_Results/',gene,'_0p2and0p4/',Sp_Type,'/Distributions/4State_k21only_noloop/Analyses.mat'];
        load(fil)
        MOD_FSP = Hog.MOD;
        switch i_nacl
            case 1
                Dat_Set = RNAdata.Step02;
                n_reps = 2;
                MOD_FSP.Salt = 0.2;
                nacl = '0p2NaCl';
                dat_col = [0.8 0 0.8];
                mod_col = [0 0 0];
                fil_SSA = ['Fit_Results/SSATS_',gene,'_Distributions_0.2_',Sp_Type,'.mat'];
            case 2
                Dat_Set = RNAdata.Step04;
                n_reps = 3;
                nacl = '0p4NaCl';
                MOD_FSP.Salt = 0.4;
                dat_col = [0 0.8 0];
                mod_col = [0 0 0.4];
                fil_SSA = ['Fit_Results/SSATS_',gene,'_Distributions_0.4_',Sp_Type,'.mat'];
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
        for irep=1:n_reps
            %% Step 1 -- Find median intensity of mature mRNA spot.
            if i_nacl==1
                N_cells=RNAFISHexp(1).rep(irep).NoCells;
            elseif i_nacl==2
                arep = [2 3 1];
                N_cells=RNAFISHexp(2).rep(arep(irep)).NoCells;
            end
            switch gene
                case 'STL1'
                    Dat_Set_B = Dat_Set.rep(irep).STL1;
                case 'CTT1'
                    Dat_Set_B = Dat_Set.rep(irep).CTT1;
            end
            all_mRNA_spots = [];
            for it=1:length(Dat_Set_B.low.tt)
                for ic = 1:N_cells(it)
                    spots = squeeze(Dat_Set_B.low.RNAintTSNuc(ic,:,it));
                    spots = sort(spots(~isnan(spots)));
                    n_spots = length(spots);
                    if n_spots>=n&&n_spots<=nmax
                        all_mRNA_spots = [all_mRNA_spots,spots(1:end-n)];
                    end
                end
            end
            figure(50)
            [H1,X] = hist(all_mRNA_spots,linspace(0,3e5,300));
            stairs(X,H1/sum(H1));hold on;
            md = median(all_mRNA_spots);
            
            %% Step 2a -- Find the Nuclear TS spots that are >=N*md (DATA)
            TS_Dat = zeros(max(N_cells),length(N_cells));
            for it=1:length(Dat_Set_B.low.tt)
                for ic = 1:N_cells(it)
                    spots = squeeze(Dat_Set_B.low.RNAintTSNuc(ic,:,it));
                    brightest = max(spots(~isnan(spots)));
                    if ~isempty(brightest)
                        TS_Dat(ic,it) = brightest/md;
                    else
                        TS_Dat(ic,it) = 0;
                    end
                end
            end
            %% Step 3a -- Plot the fraction of cells with TS versus time (DATA).
            
            N_spots = sum(TS_Dat>N);
            
            frac_TS_v_t = sum(TS_Dat>N)./N_cells';
            f2=figure(10*(i_nacl-1)+2);%subplot(3,1,1);
            plot(Dat_Set_B.low.tt,frac_TS_v_t,'o','color',dat_col,'markerfacecolor',dat_col,'markersize',14); hold on;
            set(gca,'fontsize',14,'ylim',[0 1])
            xlabel('time (min)')
            ylabel('fraction with active TS')
            
            %% Step 5a -- Plot the TS distributions versus time (DATA).
            f3=figure(10*(i_nacl-1)+3);
            for it=1:length(Dat_Set_B.low.tt)
                if length(Dat_Set_B.low.tt)==16||it==1
                    subplot(4,4,it);
                else
                    subplot(4,4,it+1);
                end
                TS_v = TS_Dat(1:N_cells(it),it);
                [A,B] = hist(TS_v,[N:40]);
                stairs(B,A/sum(A),'linewidth',2,'color',dat_col); hold on
                set(gca,'ylim',[0 0.05])
            end
            
            %% Step 6a -- Plot the TS distributions averaged over time (DATA).
            f4=figure(14+i_nacl);
            Ts_all=[];
            for it=1:length(Dat_Set_B.low.tt)
                TS_v = TS_Dat(1:N_cells(it),it);
                Ts_all = [Ts_all,TS_v'];
                TS_all_reps = [TS_all_reps,Ts_all];
                
                TS_All(i_nacl).Data = [TS_All(i_nacl).Data,TS_v'];
                
            end
            [A,B] = hist(Ts_all,[N:40]);
            stairs(B,A/sum(A),'linewidth',2); hold on
            set(gca,'fontsize',14)
            xlabel('mRNA per active TS')
            ylabel('Probability')
            
            Means_Active_TS(ig,i_nacl,irep) = mean(Ts_all(Ts_all>N));
            Median_Active_TS(ig,i_nacl,irep) = median(Ts_all(Ts_all>N));
            
            if irep==1
                if length(N_cells)==16
                    Means_Active_TS_Mod(ig,i_nacl)=Model_TS_Results([],MOD_FSP,[],PD,N,f2,f3,f4,mod_col,max(N_cells),N_cells,inten_dist,121);                    
                else
                    Means_Active_TS_Mod(ig,i_nacl)=Model_TS_Results([],MOD_FSP,[],PD,N,f2,f3,f4,mod_col,max(N_cells),N_cells([1,1:15]),inten_dist,121);
                end
            end
        end
        %% Save figures
        set(f4,'position',[1245         809         274         161]);
        figure(f4);set(gca,'fontsize',14,'yscale','log','ytick',10.^[-4:0],'ylim',[1e-3,1],'xlim',[0,30])
        xlabel('')
        ylabel('')
        
        set(f2,'position',[1245         809         274         161]);
        figure(f2);set(gca,'fontsize',14,'yscale','linear','xlim',[0,60])
        ylabel('')
        xlabel('')
 
        saveas(f2,['Figures/1_TS_Figures/',gene,'_',nacl,'_',Sp_Type,'_TS_Fraction_vs_t'],'fig')
        saveas(f4,['Figures/1_TS_Figures/',gene,'_',nacl,'_',Sp_Type,'_TS_Ave_Distributions'],'fig')
    end
    
end




