function Make_Panels_For_Figure_4(GENE)
SPACE = 'Spatial';  %'Spatial' or 'NonSpatial'
% DIST = {'Means','Moments_Corrected','Distributions'};
%% Plot MCMC Ellipses
f41=figure(41);clf;
set(f41,'position',[1000        1195         325         143]);
f1= figure(1);
close(f1);
uiopen(['Figures/MCMC_',GENE,'_Spatial_All_Pars_Reformatted.fig'],1);
f1=figure(1);

figure(f41);
copyobj(f1.Children(2).Children([1 2 7:16]),gca);
close(f1);

uiopen(['Figures/MCMC_',GENE,'_NonSpatial_All_Pars_Reformatted.fig'],1);
f1=figure(1);

figure(f41);
copyobj(f1.Children(2).Children([1 2 7:16]),gca);
close(f1);

pla = get(gca,'Children');
set(pla(2),'XData',[-4,4])
set(pla(1),'YData',[-4,5])
set(gca,'xlim',[-2.5,4],'ylim',[-4,5],'fontsize',15,'ytick',[-4 -2 0 2 4])
box on

pca = get(gca,'Children');
set(pca([4:2:12,13]),'Visible','off')
set(pca([3,15]),'color','m') % vars 4th
pca(16).FaceColor = 'm';
set(pca([5,17]),'color',[0 0 0]) % dists
pca(18).FaceColor = [0 0 0];
set(pca([7,19]),'color',[0 0 1]) % vars-xi
pca(20).FaceColor = [0 0 1];
set(pca([9,21]),'color',[1 0.0 0.0]) % means
pca(22).FaceColor = [1 0.0 0.0];

%% Parameter uncertainty
f42 = figure(42);clf
set(f42,'position',[999   858   272   141]);

uiopen(['Figures/MCMC_',GENE,'_Errors_Red_Pars.fig'],1);

f1= figure(1);
A = [f1.Children(2).Children(3).YData;...
     f1.Children(2).Children(4).YData]';
B = [f1.Children(2).Children(5).YData;...
     f1.Children(2).Children(6).YData]';

figure(f42);
bar([1 2.1 3.2 4.3],A([1 2 4 3],2),0.45,'r'); hold on;
bar([1.5 2.6 3.7 4.8],A([1 2 4 3],1),0.45,'r'); hold on;
bar([1 2.1 3.2 4.3],B([1 2 4 3],2),0.45,'b'); hold on;
bar([1.5 2.6 3.7 4.8],B([1 2 4 3],1),0.45,'y'); hold on;
set(gca,'yscale','log','fontsize',15,'ylim',[1e-6 1e15],'ytick',[1e-10 1e-5 1e0 1e5 1e10 1e15])
set(gca,'xtick',[1.25 2.35 3.45 4.55],'xticklabel',{'\mu(t)','\mu(t),\sigma(t)','4th','P(t)'})
set(gca,'xlim',[0.7 5.1])
% 
% figure(f42);
% bar([1 2.1 3.2 4.3 5.4],A([1 2 4 5 3],2),0.45,'r'); hold on;
% bar([1.5 2.6 3.7 4.8 5.9],A([1 2 4 5 3],1),0.45,'r'); hold on;
% bar([1 2.1 3.2 4.3 5.4],B([1 2 4 5 3],2),0.45,'b'); hold on;
% bar([1.5 2.6 3.7 4.8 5.9],B([1 2 4 5 3],1),0.45,'y'); hold on;
% set(gca,'yscale','log','fontsize',15,'ylim',[1e-6 1e15],'ytick',[1e-10 1e-5 1e0 1e5 1e10 1e15])
% set(gca,'xtick',[1.25 2.35 3.45 4.55 5.65],'xticklabel',{'\mu(t)','\mu(t),\sigma(t)','4th','Gaussian','P(t)'})
% set(gca,'xlim',[0.7 6.2])
% 
% Parameter bias
f43 = figure(43);clf
set(f43,'position',[999   658   272   141]);

f1= figure(1);
A = [f1.Children(4).Children(3).YData;...
    f1.Children(4).Children(4).YData]';
B = [f1.Children(4).Children(5).YData;...
    f1.Children(4).Children(6).YData]';

figure(f43);
bar([1 2.1 3.2 4.3],A([1 2 4 3],2),0.45,'r'); hold on;
bar([1.5 2.6 3.7 4.8],A([1 2 4 3],1),0.45,'r'); hold on;
bar([1 2.1 3.2 4.3],B([1 2 4 3],2),0.45,'b'); hold on;
bar([1.5 2.6 3.7 4.8],B([1 2 4 3],1),0.45,'y'); hold on;
set(gca,'yscale','log','fontsize',15,'ylim',[1e-6 1e15],'ytick',[1e-10 1e-5 1e0 1e5 1e10 1e15])
set(gca,'xtick',[1.25 2.35 3.45 4.55],'xticklabel',{'\mu(t)','\mu(t),\sigma(t)','4th','P(t)'})
set(gca,'xlim',[0.7 5.1])

% figure(f43);
% bar([1 2.1 3.2 4.3 5.4],A([1 2 4 5 3],2),0.45,'r'); hold on;
% bar([1.5 2.6 3.7 4.8 5.9],A([1 2 4 5 3],1),0.45,'r'); hold on;
% bar([1 2.1 3.2 4.3 5.4],B([1 2 4 5 3],2),0.45,'b'); hold on;
% bar([1.5 2.6 3.7 4.8 5.9],B([1 2 4 5 3],1),0.45,'y'); hold on;
% set(gca,'yscale','log','fontsize',15,'ylim',[1e-6 1e15],'ytick',[1e-10 1e-5 1e0 1e5 1e10 1e15])
% set(gca,'xtick',[1.25 2.35 3.45 4.55 5.65],'xticklabel',{'\mu(t)','\mu(t),\sigma(t)','4th','Gaussian','P(t)'})
% set(gca,'xlim',[0.7 6.2])

close(f1);

%% TS Average Loading 
f44 = figure(44);clf
set(f44,'position',[1025         773         273         155]);
uiopen(['Figures/1_TS_Figures/Mean_Intensity_Different_Analyses.fig'],1);
f1=figure(1);
if strcmp(GENE,'STL1')     %'STL1' or 'CTT1'
    figure(f1);
    subplot(2,1,1);
else
    figure(f1);
    subplot(2,1,2);
end
ax = gca;

figure(f44);
copyobj(allchild(ax),gca);
set(gca,'yscale','log','fontsize',14,'ylim',[1e0 1e8],'ytick',10.^[0:2:8],'xlim',[0.5 5.3])
set(gca,'xtick',[1.25 2.35 3.45 4.55],'xticklabel',{'\mu(t)','\mu(t),\sigma(t)','4th','P(t)'})
box on
close(1)

%% Plot TS fraction vs time.
nacl = '0p2NaCl';
f45 = figure(45);clf
set(f45,'position',[725   565   261   163]);
uiopen(['Figures/1_TS_Figures/',GENE,'_',nacl,'_',SPACE,'_TS_Fraction_vs_t.fig'],1);
figure(1);
ax = gca;
figure(f45);
copyobj(allchild(ax),gca);
figure(f45);set(gca,'fontsize',14,'yscale','linear','xlim',[0,60])
box on
close(1)

nacl = '0p4NaCl';
f46 = figure(46);clf
set(f46,'position',[725   365   261   163]);
uiopen(['Figures/1_TS_Figures/',GENE,'_',nacl,'_',SPACE,'_TS_Fraction_vs_t.fig'],1);
figure(1);
ax = gca;
figure(f46);
copyobj(allchild(ax),gca);
figure(f46);set(gca,'fontsize',14,'yscale','linear','xlim',[0,60])
box on
close(1)

%% Plot TS Distributions
nacl = '0p2NaCl';
f47 = figure(47);clf
set(f47,'position',[1025         573         273         155]);
uiopen(['Figures/1_TS_Figures/',GENE,'_',nacl,'_',SPACE,'_TS_Ave_Distributions.fig'],1);
figure(1);
ax = gca;
figure(f47);
copyobj(allchild(ax),gca);
figure(f47);set(gca,'fontsize',14,'yscale','log','ytick',10.^[-4:0],'ylim',[1e-3,1],'xlim',[0,30])
box on
close(1)

nacl = '0p4NaCl';
f48 = figure(48);clf
set(f48,'position',[1025         373         273         155]);
uiopen(['Figures/1_TS_Figures/',GENE,'_',nacl,'_',SPACE,'_TS_Ave_Distributions.fig'],1);
figure(1);
ax = gca;
figure(f48);
copyobj(allchild(ax),gca);
figure(f48);set(gca,'fontsize',14,'yscale','log','ytick',10.^[-4:0],'ylim',[1e-3,1],'xlim',[0,30])
box on
close(1)

%%
if strcmp(GENE,'CTT1')
    figure(41);
    set(gca,'ylim',[-12 4],'xlim',[-3.5 1],'ytick',(-12:4:4));
    pca = get(gca,'Children');
    pca(1).YData=[-12,4];
    pca(2).XData=[-4,4];
    
    figure(42);
    set(gca,'ylim',[1e-6 1e15],'ytick',10.^(-5:5:25))
    
    figure(43);
    set(gca,'ylim',[1e-6 1e22],'ytick',10.^(-0:10:25))
    
    figure(44);
    set(gca,'ylim',[1e-2 1e8],'ytick',10.^[-2:2:8])
    
    figure(45);
    set(gca,'ylim',[0 0.8],'ytick',[0 0.2 0.4 0.6 0.8])
    
    figure(46);
    set(gca,'ylim',[0 0.8],'ytick',[0 0.2 0.4 0.6 0.8])
       
end

