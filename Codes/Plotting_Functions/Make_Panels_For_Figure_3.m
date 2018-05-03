function Make_Panels_For_Figure_3(gene,nacl)
%% SSA Figures
% gene = 'STL1';
% nacl = 0.2;
% fil = ['Fit_Results/SSA_',gene,'_',num2str(nacl),'.mat'];
fil =['Fit_Results/',gene,'_0p2and0p4/NonSpatial/Distributions/4State_k21only_noloop/Analyses.mat'];
load(fil,'Hog')
switch(nacl)
    case 0.2
        SSA = Hog.SSA(1);
        MOD = Hog.MOD;
        MOD.Salt = 0.2;
    case 0.4
        SSA = Hog.SSA(2);
        MOD = Hog.MOD;
        MOD.Salt = 0.4;
end
MOD.SSA_Expts = 100; % Number of SSA experiments.
MOD.SSA_Cells = 200; % Number of cells per experiment.

%%  Generate data and make rough plots.
MOD.f = 0.00;  % Probability of spot ovrlap.
MOD.plot_SSA_dist_results(MOD,SSA,32,122,1,[])
MOD.plot_SSA_moment_results(MOD,SSA,124,125)
MOD.Type='Moments';
MOD.Fit_Type='Moments';
MOD.Moment_Order = 2;
Moments = Hog_Model.solve_hog(MOD);
MOD.plot_moments_mod(MOD,Moments,124,125)
MOD.Type='Distributions';
Dists = Hog_Model.solve_hog(MOD);
MOD.plot_distributions_mod(MOD,Dists,32,122)
for TMP = {'_A','_B','_C'}
    MOD.Type = 'Distributions';
    MOD.Reps = TMP{1};
    MOD.plot_distributions_dat(MOD,MOD.Data,32,122)   
    MOD.Type = 'Moments';
    MOD.plot_moments_dat(MOD,MOD.Data,124,125)
end

%% Reformat distribution plots
f32=figure(32);
COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4; 0.8 0.8 0.8 ];
COLa = ['r';'b';'c'];
set(f32,'position',[1928 290  409 385])
for j=1:16
    subplot(4,4,j)
    set(gca,'xlim',[0 150],'ylim',[0 0.04],'yticklabel',[],'xticklabel',[],'fontsize',13);
    
    pla = get(gca,'children');
    if nacl==0.2
        for i=1:3
            set(pla(i),'color',COLa(i),'linewidth',2);
        end
        set(pla(5),'color','m','linewidth',3)
        set(pla(4),'color','k','linewidth',3);
        for i=6:length(pla)
            set(pla(i),'color',[0.6 0.6 0.6],'linewidth',1);
        end
    else
        for i=1:3
            set(pla(i),'color',COLa(i),'linewidth',2);
        end
        set(pla(5),'color','m','linewidth',3)
        set(pla(4),'color','k','linewidth',3);
        for i=6:length(pla)
            set(pla(i),'color',[0.6 0.6 0.6],'linewidth',1);
        end
    end
end
for j=1:16
    %             text(70,0.035,MOD.Gene,'fontsize',14);
    %             text(70,0.03,[num2str(nacl),'M NaCl'],'fontsize',14);
    subplot(4,4,j)
    text(70,0.035,[num2str(MOD.tt(j)),' min'],'fontsize',14);
    x = mod(j-1,4);y = floor((j-1)/4);
    set(gca,'Position',[0.01+x/4 0.77-y/4 0.23 0.23])
end

%% update SSA summary stat plots
f31 = figure(31);clf;
set(f31,'position',[891   432   290   442]);
figure(124); ax = gca;

figure(31);subplot(3,1,1); hold off
copyobj(allchild(ax),gca);
set(gca,'fontsize',15,'xticklabel',[],'ylim',[0.01 100],'ytick',[0.01 0.1 1 10 100],...
    'yscale','log','yticklabel',{'0.01','0.1','1','10','100'})
pla = get(gca,'children');ylabel('\mu_{mRNA}')
set(pla(1),'markersize',14,'marker','^','linestyle','none','markerfacecolor','b')
set(pla(2),'markersize',14,'marker','^','linestyle','none','markerfacecolor','r')
set(pla(3),'markersize',14,'marker','^','linestyle','none','markerfacecolor','c')
set(pla(4),'color','k','linewidth',3)
set(pla(5),'color','m','linewidth',3)
set(pla(6:end),'color',[0.6 0.6 0.6])
pla(1).Visible='off';
%         title('Summary Statistics')
xlabel('');
ylabel('')
set(gca,'xlim',[0 55],'xlim',[0 55],'xtick',[0:10:50],'xticklabel',[0:10:50])

figure(125); ax = gca;
figure(31);subplot(3,1,2); hold off
copyobj(allchild(ax),gca);
set(gca,'fontsize',15,'xticklabel',[],'ylim',[0.1 100],'ytick',[0.1 1 10 100],...
    'yscale','log','yticklabel',{'0.1','1','10','100'})
pla = get(gca,'children');ylabel('\sigma_{mRNA}')
set(pla(1),'markersize',14,'marker','^','linestyle','none','markerfacecolor','b')
set(pla(2),'markersize',14,'marker','^','linestyle','none','markerfacecolor','r')
set(pla(3),'markersize',14,'marker','^','linestyle','none','markerfacecolor','c')
set(pla(4),'color','k','linewidth',3)
set(pla(5),'color','m','linewidth',3)
set(pla(6:end),'color',[0.6 0.6 0.6])
pla(1).Visible='off';
xlabel('');
ylabel('')
set(gca,'xlim',[0 55],'xlim',[0 55],'xtick',[0:10:50],'xticklabel',[0:10:50])

figure(122); subplot(1,2,1); ax = gca;
figure(31);subplot(3,1,3); hold off
copyobj(allchild(ax),gca);
set(gca,'fontsize',15,'xticklabel',[],'ylim',[0.001 1],'ytick',[0.001 0.01 0.1 1 ],...
    'yscale','log','yticklabel',{'0.001','0.01','0.1','1'})
pla = get(gca,'children');ylabel('P(ON)')
set(pla(1),'markersize',14,'marker','^','linestyle','none','markerfacecolor','r')
set(pla(2),'markersize',14,'marker','^','linestyle','none','markerfacecolor','b')
set(pla(3),'markersize',14,'marker','^','linestyle','none','markerfacecolor','c')
set(pla(4),'color','k','linewidth',3)
set(pla(5),'color','m','linewidth',3)
set(pla(6:end),'color',[0.6 0.6 0.6])
pla(2).Visible='off';
xlabel('');
ylabel('')
set(gca,'xlim',[0 55],'xlim',[0 55],'xtick',[0:10:50],'xticklabel',[0:10:50],'ylim',[0.005 1])

close(122);close(125);close(124);