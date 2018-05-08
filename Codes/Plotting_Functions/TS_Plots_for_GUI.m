% To find the "standard" mRNa intensity, I look at all cells at all time
% points. I will then rank them from the dimmest to the brightest.  Then I
% will remove the 'n' brightest spots (in case these are TS).  Finally, I
% will find the median value of the remaining spots.
% clear all
speed = GUI.mRNAelongationEditField.Value; % CTT1  - mn of on TS.
Sp_Type = GUI.MOD.Spatial;

if isempty(GUI.RNAdata)
    try
        disp('Loading Data -- This could take a few minutes the first time.');
        load('Hog_Data/RNAspots.mat')
        disp('Loading data complete');
    catch
        msgbox({'Missing data files that are too large for GitHub.',...
            'Download contennts of the folder:',...
            'https://www.dropbox.com/sh/88sbqt5g5qfe5o5/AAAEnqOem_t5E1mMoshXHnv8a?dl=0',...
            'and place within the Hog_Data folder in the main directory.'})
    end
        GUI.RNAdata = RNAdata;
    clear RNAdata
end
load 'Hog_Data/Analyzed_FISH_Data_2016_06';
load('Hog_Data/Hog_Probes.mat')

n=3;
nmax = inf;

fignum = GUI.FigEditField_4.Value;
% Number of mRNA spots required to positively ID as an mRNA. In this
% analysis, I am only going to look at cells that have more than 'n' spots
% and then use these to determine what the average mature mRNA looks like.
% This is meant to avoid using the TS when trying to determine what is a
% mature mRNA.
N = 2; % Minimum number of nascent mRNA required to be considered as a TS.
clear TMP TS_All
switch GUI.GeneTypeButtonGroup.SelectedObject.Text
    case 'STL1'
        gene = 'STL1';
        PD = PD_STL1;
        len = 1630;
    case 'CTT1'
        gene = 'CTT1';
        PD = PD_CTT1;
        len  = 1618;
end

for i_nacl = 1:2
    TS_All(i_nacl).Data = [];
    MOD_FSP = GUI.MOD;
    switch i_nacl
        case 1
            Dat_Set = GUI.RNAdata.Step02;
            n_reps = 2;
            MOD_FSP.Salt = 0.2;
            nacl = '0p2NaCl';
            dat_col = [0.8 0 0.8;...
                       0 0.8 0;...
                       0 0 0.8];
            mod_col = [0 0 0];
        case 2
            Dat_Set = GUI.RNAdata.Step04;
            n_reps = 3;
            nacl = '0p4NaCl';
            MOD_FSP.Salt = 0.4;
            dat_col = [0.8 0 0.8;...
                0 0.8 0;...
                0 0 0.8];            
            mod_col = [0 0 0.4];
    end
    All_TS_Msmts = [];
    
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
%         figure(fignum+i_nacl)
%         [H1,X] = hist(all_mRNA_spots,linspace(0,3e5,300));
%         stairs(X,H1/sum(H1));hold on;
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
        f2=figure(fignum+10*(i_nacl)+2);%subplot(3,1,1);
        plot(Dat_Set_B.low.tt,frac_TS_v_t,'o','color',dat_col(irep,:),'markerfacecolor',dat_col(irep,:),'markersize',14); hold on;
        set(gca,'fontsize',14,'ylim',[0 1])
        xlabel('time (min)')
        ylabel('fraction with active TS')
      
        %% Step 5a -- Plot the TS distributions versus time (DATA).
        f3=figure(fignum+10*(i_nacl)+3);
        for it=1:length(Dat_Set_B.low.tt)
            if length(Dat_Set_B.low.tt)==16||it==1
                subplot(4,4,it);
            else
                subplot(4,4,it+1);
            end
            TS_v = TS_Dat(1:N_cells(it),it);
            [A,B] = hist(TS_v,[N:40]);
            stairs(B,A/sum(A),'s','linewidth',2,'color',dat_col(irep,:)); hold on
            plot([0 30],1/N_cells(it)*[1,1],'--','color',dat_col(irep,:));
        end
        
        fnas = figure(fignum+37+i_nacl);
        set(fnas,'position',[1561     1218-200*i_nacl         560         100]);
        figure(fnas);
        tar = [5 6 7 9 11];
        for it=1:length(tar)
            subplot(1,length(tar),it);
            if length(Dat_Set_B.low.tt)==16||it==1
                itt = tar(it);
            else
                itt = tar(it)-1;
            end
            TS_v = TS_Dat(1:N_cells(itt),itt);
            [A,B] = hist(TS_v,[N:40]);
            stairs(B,A/sum(A),'s','linewidth',2,'color',dat_col(irep,:)); hold on
            plot([0 30],1/N_cells(itt)*[1,1],'--','color',dat_col(irep,:));
            set(gca,'yscale','log')
        end

        
        %% Step 6a -- Plot the TS distributions averaged over time (DATA).
        f4=figure(fignum+2+i_nacl);
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
        
        Means_Active_TS(i_nacl,irep) = mean(Ts_all(Ts_all>N));
        Median_Active_TS(i_nacl,irep) = median(Ts_all(Ts_all>N));
        All_TS_Msmts = [All_TS_Msmts;A];
        %% Step 2b -- Find the TS spots intensities (MODEL)

            if length(N_cells)==16
                Means_Active_TS_Mod(i_nacl)=Model_TS_Results([],MOD_FSP,[],PD,N,f2,f3,f4,mod_col,max(N_cells),N_cells,inten_dist,fnas);
            else
                Means_Active_TS_Mod(i_nacl)=Model_TS_Results([],MOD_FSP,[],PD,N,f2,f3,f4,mod_col,max(N_cells),N_cells([1,1:15]),inten_dist,fnas);
            end
    end
        A = sum(All_TS_Msmts);
        stairs(B,A/sum(A),'linewidth',2); hold on
        set(gca,'fontsize',14)
        xlabel('mRNA per active TS')
        ylabel('Probability')

    %% Save figures
    set(f4,'position',[145         809         274         161]);
    figure(f4);set(gca,'fontsize',14,'yscale','log','ytick',10.^[-4:0],'ylim',[1e-3,1],'xlim',[0,30])
    xlabel('')
    ylabel('')
    
    set(f2,'position',[145         809         274         161]);
    figure(f2);set(gca,'fontsize',14,'yscale','linear','xlim',[0,60])
    ylabel('')
    xlabel('')
    
end





