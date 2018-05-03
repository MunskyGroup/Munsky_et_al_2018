%% Run codes to load data and compute the full FSP model TS predictions.
if ~exist('Mn_Mod_NSP','var')
    Sp_Type = 'NonSpatial';
    Find_TS_in_Data;
    Mn_Mod_NSP = mean(Means_Active_TS_Mod,2);
end

if ~exist('Mn_Mod_SP','var')
    Sp_Type = 'Spatial';
    Find_TS_in_Data;
    Mn_Mod_SP = mean(Means_Active_TS_Mod,2);
end

%% Pull parameters from mature mRNA fits
f44 = figure(44);
set(f44,'position',[1000         886         286         366]);
A = readtable('FOUND_PARS_Red_Pars.csv');

clear Kr Kv
Kr(:,1) = A.Var27;
Kr(:,2) = A.Var31;
Kr(:,3) = A.Var35;
Kv(:,1) = A.Var29;
Kv(:,2) = A.Var33;
Kv(:,3) = A.Var37;

L_CTT1 = 1618; %(length of region with CTT1 probes)
L_STL1 = 1630; %(length of region with STL1 probes)

%%  Use simplified theory to compute average elongation rates.
TS_Ctt1 = 1/5*sum(sum(squeeze(Means_Active_TS(2,:,:))));
STD_CTT1 = sqrt(1/5*sum(sum(squeeze(Means_Active_TS(2,:,:).^2))) -TS_Ctt1^2);
ke_CTT1 = 1/2./TS_Ctt1*L_CTT1*Kr(8,2);
ke_CTT1_arr = 1/2./squeeze(Means_Active_TS(2,:,:))*L_CTT1*Kr(8,2);

TS_STL1 = 1/5*sum(sum(squeeze(Means_Active_TS(1,:,:))));
STD_STL1 = sqrt(1/5*sum(sum(squeeze(Means_Active_TS(1,:,:).^2))) -TS_STL1^2);
ke_STL1 = 1/2./TS_STL1*L_STL1*Kr(18,2);
ke_STL1_arr = 1/2./squeeze(Means_Active_TS(1,:,:))*L_STL1*Kr(18,2);

%% Find CTT1 Elongation rates -- Simplified Theory.
for i = 1:10
    AA = zeros(5,1);
    tmp = 1/2./squeeze(Means_Active_TS(1,:,:))*L_CTT1*max(Kr(i,1:3),[],2);
    AA(:) = tmp(isfinite(tmp));
    mn(i,1) = mean(AA);
    sem(i,1) = std(AA)/sqrt(5);
end

%% Find STL1 Elongation rates -- Simplified Theory.
for i = 1:10
    AA = zeros(5,1);
    tmp = 1/2./squeeze(Means_Active_TS(2,:,:))*L_STL1*max(Kr(i+10,1:3),[],2);
    AA(:) = tmp(isfinite(tmp));
    mn(i,1) = mean(AA);
    sem(i,1) = std(AA)/sqrt(5);
end

%% Make plots of mRNA/TS for simplified theory (first 8 bars) and full simulations (last 2 bars).
ke = 1/2/TS_Ctt1*L_CTT1*Kr(8,2);
TS(1:10,:) = 1/2/ke*L_CTT1.*Kr(1:10,:);
TS(11:20,:) = 1/2/ke*L_STL1.*Kr(11:20,:);

AA1 = [max(TS(1,:)),max(TS(6,:));... % means
       max(TS(2,:)),max(TS(7,:));... % means and vars
       max(TS(4,:)),max(TS(9,:));... % extended
       max(TS(3,:)),max(TS(8,:));... % FSP
       Mn_Mod_NSP(2),Mn_Mod_SP(2)];

subplot(2,1,2);
bar([1 2.1 3.2 4.3],AA1([1 2 3 5],1),0.43,'b'); hold on;
bar([1.5 2.6 3.7 4.8],AA1([1 2 3 5],2),0.43,'y'); hold on;
set(gca,'yscale','log','fontsize',14,'ylim',[1e-2 1e8],'ytick',10.^[-2:2:8],'xlim',[0.5 5.3])
set(gca,'xtick',[1.25 2.35 3.45 4.55],'xticklabel',{'\mu(t)','\mu(t),\sigma(t)','4th','P(t)'})
% title('CTT1 Average mRNA per active TS')
hold on
plot([0.5 5.5],Means_Active_TS(2,1,1)*[1 1],'m--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(2,1,2)*[1 1],'m--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(2,2,1)*[1 1],'c--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(2,2,2)*[1 1],'c--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(2,2,3)*[1 1],'c--','linewidth',2)


AA2 = [max(TS(11,:)),max(TS(16,:));...
       max(TS(12,:)),max(TS(17,:));...
       max(TS(14,:)),max(TS(19,:));...
       max(TS(13,:)),max(TS(18,:));...
       Mn_Mod_NSP(1),Mn_Mod_SP(1)];

subplot(2,1,1);
bar([1 2.1 3.2 4.3],AA2([1 2 3 5],1),0.43,'b'); hold on;
bar([1.5 2.6 3.7 4.8],AA2([1 2 3 5],2),0.43,'y'); hold on;
set(gca,'yscale','log','fontsize',14,'ylim',[1e-2 1e8],'ytick',10.^[-2:2:8],'xlim',[0.5 5.3])
set(gca,'xtick',[1.25 2.35 3.45 4.55],'xticklabel',{'\mu(t)','\mu(t),\sigma(t)','4th','P(t)'})
% title('STL1 Average mRNA per active TS')
hold on
plot([0.5 5.5],Means_Active_TS(1,1,1)*[1 1],'m--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(1,1,2)*[1 1],'m--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(1,2,1)*[1 1],'c--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(1,2,2)*[1 1],'c--','linewidth',2)
plot([0.5 5.5],Means_Active_TS(1,2,3)*[1 1],'c--','linewidth',2)


saveas(f44,['Figures/1_TS_Figures/Mean_Intensity_Different_Analyses'],'fig')


%% 
disp('********TS analysis results*********')

mn_STL1_TS_exp = 1/5*sum(sum(squeeze(Means_Active_TS(1,:,:))))
range_STL1_TS_exp(1:6) = Means_Active_TS(1,:,:)
Mod_Pred_TS_mn = AA2(1,:)
Mod_Pred_TS_mn_var = AA2(2,:)
Mod_Pred_TS_ext = AA2(3,:)
Mod_Pred_TS_FSP = AA2(5,:)


mn_CTT1_TS_exp = 1/5*sum(sum(squeeze(Means_Active_TS(2,:,:))))
range_CTT1_TS_exp(1:6) = Means_Active_TS(2,:,:)
Mod_Pred_TS_mn = AA1(1,:)
Mod_Pred_TS_mn_var = AA1(2,:)
Mod_Pred_TS_ext = AA1(3,:)
Mod_Pred_TS_FSP = AA1(5,:)

