%% Update_All_Figures
clear all
close all
path(path,'Plotting_Functions');
path(path,'Hog_Data');
recompute = 0;
re_mcmc = 0;
refig = 0;
if recompute==1 % Recompute all fit results
    %% Delete all old figures
    delete Figures/*.fig
    delete Figures/0_Main_Figures/*.fig
    delete Figures/1_TS_Figures/*.fig
    delete Final_Figures/*/*.fig
    delete Final_Figures/*/*.jpg
    refig=1;re_mcmc=1;    
    %% Update all of the analyses.
    Update_Analyses  %-- Only run after parameters have been changed.
end

if re_mcmc==1 % Recompile all MCMC results
    Update_MCMC      %-- Only run after parameters have been changed.
    refig=1;
end

    %% Prepare draft figures
if refig==1 % Remake draft figures.
    Update_MCMC_Plots
    Update_Plots
    Plot_TS_Theory
end
%% Reformat to Figure 2 
close all
Make_Panels_For_Figure_2('STL1')
saveas(21,'Final_Figures/Fig2/2A','fig');saveas(21,'Final_Figures/Fig2/2A','epsc')
saveas(22,'Final_Figures/Fig2/2B','fig');saveas(22,'Final_Figures/Fig2/2B','epsc')
saveas(23,'Final_Figures/Fig2/2C','fig');saveas(23,'Final_Figures/Fig2/2C','epsc')
saveas(24,'Final_Figures/Fig2/2D','fig');saveas(24,'Final_Figures/Fig2/2D','epsc')
saveas(121,'Final_Figures/SupMat/S10/S10AB','fig');saveas(121,'Final_Figures/SupMat/S10/S10AB','epsc')
saveas(122,'Final_Figures/SupMat/S10/S10CD','fig');saveas(122,'Final_Figures/SupMat/S10/S10CD','epsc')

%% Reformat to Figure 3
Make_Panels_For_Figure_3('STL1',0.2)
saveas(31,'Final_Figures/Fig3/3ABC','fig');saveas(31,'Final_Figures/Fig3/3ABC','epsc')
saveas(32,'Final_Figures/Fig3/3D','fig');saveas(32,'Final_Figures/Fig3/3D','epsc')

%% Reformat to Figure 4
Make_Panels_For_Figure_4('STL1')
saveas(41,'Final_Figures/Fig4/4A','fig');saveas(41,'Final_Figures/Fig4/4A','epsc')
saveas(42,'Final_Figures/Fig4/4B','fig');saveas(42,'Final_Figures/Fig4/4B','epsc')
saveas(43,'Final_Figures/Fig4/4C','fig');saveas(43,'Final_Figures/Fig4/4C','epsc')
saveas(44,'Final_Figures/Fig4/4D','fig');saveas(44,'Final_Figures/Fig4/4D','epsc')
saveas(45,'Final_Figures/Fig4/4E','fig');saveas(45,'Final_Figures/Fig4/4E','epsc')
saveas(46,'Final_Figures/Fig4/4F','fig');saveas(46,'Final_Figures/Fig4/4F','epsc')
saveas(47,'Final_Figures/Fig4/4G','fig');saveas(47,'Final_Figures/Fig4/4G','epsc')
saveas(48,'Final_Figures/Fig4/4H','fig');saveas(48,'Final_Figures/Fig4/4H','epsc')

%% Reformat to Finalize Figures (Main Text)
%% S5  Hog1p vs time
path(path,'Hog_Data/Hog1p_Data_Fitting')
figure(51);Main_ID_Match_Hog1
saveas(51,'Final_Figures/SupMat/S5/S5','fig');saveas(51,'Final_Figures/SupMat/S5/S5','epsc')

%%
close all
Make_Panels_For_Figure_2('CTT1')
saveas(121,'Final_Figures/SupMat/S10/S10EF','fig');saveas(121,'Final_Figures/SupMat/S10/S10EF','epsc')
saveas(122,'Final_Figures/SupMat/S10/S10GH','fig');saveas(122,'Final_Figures/SupMat/S10/S10GH','epsc')

%% S6 Distributions for STL1 -- all analyses
close all
uiopen('Figures/STL1_NonSpatial_Fit_Means_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/STL1_NonSpatial_Fit_Means_Plot_Marg_Dists_0p4.fig',1)
uiopen('Figures/STL1_NonSpatial_Fit_Moments_Corrected_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/STL1_NonSpatial_Fit_Moments_Corrected_Plot_Marg_Dists_0p4.fig',1)
uiopen('Figures/STL1_NonSpatial_Fit_Distributions_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/STL1_NonSpatial_Fit_Distributions_Plot_Marg_Dists_0p4.fig',1)
uiopen('Figures/STL1_NonSpatial_Fit_Moments_4_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/STL1_NonSpatial_Fit_Moments_4_Plot_Marg_Dists_0p4.fig',1)
saveas(1,'Final_Figures/SupMat/S6/S6A','fig');saveas(1,'Final_Figures/SupMat/S6/S6A','epsc')
saveas(2,'Final_Figures/SupMat/S6/S6B','fig');saveas(2,'Final_Figures/SupMat/S6/S6B','epsc')
saveas(3,'Final_Figures/SupMat/S6/S6C','fig');saveas(3,'Final_Figures/SupMat/S6/S6C','epsc')
saveas(4,'Final_Figures/SupMat/S6/S6D','fig');saveas(4,'Final_Figures/SupMat/S6/S6D','epsc')
saveas(5,'Final_Figures/SupMat/S6/S6E','fig');saveas(5,'Final_Figures/SupMat/S6/S6E','epsc')
saveas(6,'Final_Figures/SupMat/S6/S6F','fig');saveas(6,'Final_Figures/SupMat/S6/S6F','epsc')
saveas(7,'Final_Figures/4th_Moments/S6E','fig');saveas(7,'Final_Figures/4th_Moments/S6E','epsc')
saveas(8,'Final_Figures/4th_Moments/S6F','fig');saveas(8,'Final_Figures/4th_Moments/S6F','epsc')

%% S7 Distributions for CTT1 -- all analyses
close all
uiopen('Figures/CTT1_NonSpatial_Fit_Means_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/CTT1_NonSpatial_Fit_Means_Plot_Marg_Dists_0p4.fig',1)
uiopen('Figures/CTT1_NonSpatial_Fit_Moments_Corrected_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/CTT1_NonSpatial_Fit_Moments_Corrected_Plot_Marg_Dists_0p4.fig',1)
uiopen('Figures/CTT1_NonSpatial_Fit_Distributions_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/CTT1_NonSpatial_Fit_Distributions_Plot_Marg_Dists_0p4.fig',1)
uiopen('Figures/CTT1_NonSpatial_Fit_Moments_4_Plot_Marg_Dists_0p2.fig',1)
uiopen('Figures/CTT1_NonSpatial_Fit_Moments_4_Plot_Marg_Dists_0p4.fig',1)
saveas(1,'Final_Figures/SupMat/S7/S7A','fig');saveas(1,'Final_Figures/SupMat/S7/S7A','epsc')
saveas(2,'Final_Figures/SupMat/S7/S7B','fig');saveas(2,'Final_Figures/SupMat/S7/S7B','epsc')
saveas(3,'Final_Figures/SupMat/S7/S7C','fig');saveas(3,'Final_Figures/SupMat/S7/S7C','epsc')
saveas(4,'Final_Figures/SupMat/S7/S7D','fig');saveas(4,'Final_Figures/SupMat/S7/S7D','epsc')
saveas(5,'Final_Figures/SupMat/S7/S7G','fig');saveas(5,'Final_Figures/SupMat/S7/S7G','epsc')
saveas(6,'Final_Figures/SupMat/S7/S7H','fig');saveas(6,'Final_Figures/SupMat/S7/S7H','epsc')
saveas(7,'Final_Figures/SupMat/S7/S7E','fig');saveas(7,'Final_Figures/SupMat/S7/S7E','epsc')
saveas(8,'Final_Figures/SupMat/S7/S7F','fig');saveas(8,'Final_Figures/SupMat/S7/S7F','epsc')

%% S8 Spatial Distributions for STL1 -- all analyses
close all
uiopen('Figures/STL1_Spat_Dists_Dat_0p2.fig',1)
uiopen('Figures/STL1_Spat_Dists_Dat_0p4.fig',1)
uiopen('Figures/STL1_Spatial_Fit_Distributions_Plot_Spat_Dists_Mod_0p2.fig',1)
uiopen('Figures/STL1_Spatial_Fit_Distributions_Plot_Spat_Dists_Mod_0p4.fig',1)
uiopen('Figures/STL1_Spatial_Fit_Moments_Corrected_Plot_Spat_Dists_Mod_0p2.fig',1)
uiopen('Figures/STL1_Spatial_Fit_Moments_Corrected_Plot_Spat_Dists_Mod_0p4.fig',1)
uiopen('Figures/STL1_Spatial_Fit_Means_Plot_Spat_Dists_Mod_0p2.fig',1) 
uiopen('Figures/STL1_Spatial_Fit_Means_Plot_Spat_Dists_Mod_0p4.fig',1)
uiopen('Figures/STL1_Spatial_Fit_Moments_4_Plot_Spat_Dists_Mod_0p2.fig',1)
uiopen('Figures/STL1_Spatial_Fit_Moments_4_Plot_Spat_Dists_Mod_0p4.fig',1)
saveas(1,'Final_Figures/SupMat/S8/S8A','fig');saveas(1,'Final_Figures/SupMat/S8/S8A','epsc')
saveas(2,'Final_Figures/SupMat/S8/S8B','fig');saveas(2,'Final_Figures/SupMat/S8/S8B','epsc')
saveas(3,'Final_Figures/SupMat/S8/S8C','fig');saveas(3,'Final_Figures/SupMat/S8/S8C','epsc')
saveas(4,'Final_Figures/SupMat/S8/S8D','fig');saveas(4,'Final_Figures/SupMat/S8/S8D','epsc')
saveas(5,'Final_Figures/SupMat/S8/S8E','fig');saveas(5,'Final_Figures/SupMat/S8/S8E','epsc')
saveas(6,'Final_Figures/SupMat/S8/S8F','fig');saveas(6,'Final_Figures/SupMat/S8/S8F','epsc')
saveas(7,'Final_Figures/SupMat/S8/S8G','fig');saveas(7,'Final_Figures/SupMat/S8/S8G','epsc')
saveas(8,'Final_Figures/SupMat/S8/S8H','fig');saveas(8,'Final_Figures/SupMat/S8/S8H','epsc')
saveas(9,'Final_Figures/SupMat/S8/S8I-Extra','fig');saveas(9,'Final_Figures/SupMat/S8/S8I-Extra','epsc')
saveas(10,'Final_Figures/SupMat/S8/S8J-Extra','fig');saveas(10,'Final_Figures/SupMat/S8/S8J-Extra','epsc')

%% S9 Spatial Distributions for CTT1 -- all analyses
close all
uiopen('Figures/CTT1_Spat_Dists_Dat_0p2.fig',1)
uiopen('Figures/CTT1_Spat_Dists_Dat_0p4.fig',1)
uiopen('Figures/CTT1_Spatial_Fit_Distributions_Plot_Spat_Dists_Mod_0p2.fig',1)
uiopen('Figures/CTT1_Spatial_Fit_Distributions_Plot_Spat_Dists_Mod_0p4.fig',1)
uiopen('Figures/CTT1_Spatial_Fit_Moments_Corrected_Plot_Spat_Dists_Mod_0p2.fig',1)
uiopen('Figures/CTT1_Spatial_Fit_Moments_Corrected_Plot_Spat_Dists_Mod_0p4.fig',1)
uiopen('Figures/CTT1_Spatial_Fit_Means_Plot_Spat_Dists_Mod_0p2.fig',1) 
uiopen('Figures/CTT1_Spatial_Fit_Means_Plot_Spat_Dists_Mod_0p4.fig',1)
uiopen('Figures/CTT1_Spatial_Fit_Moments_4_Plot_Spat_Dists_Mod_0p2.fig',1)
uiopen('Figures/CTT1_Spatial_Fit_Moments_4_Plot_Spat_Dists_Mod_0p4.fig',1)
saveas(1,'Final_Figures/SupMat/S9/S9A','fig');saveas(1,'Final_Figures/SupMat/S9/S9A','epsc')
saveas(2,'Final_Figures/SupMat/S9/S9B','fig');saveas(2,'Final_Figures/SupMat/S9/S9B','epsc')
saveas(3,'Final_Figures/SupMat/S9/S9C','fig');saveas(3,'Final_Figures/SupMat/S9/S9C','epsc')
saveas(4,'Final_Figures/SupMat/S9/S9D','fig');saveas(4,'Final_Figures/SupMat/S9/S9D','epsc')
saveas(5,'Final_Figures/SupMat/S9/S9E','fig');saveas(5,'Final_Figures/SupMat/S9/S9E','epsc')
saveas(6,'Final_Figures/SupMat/S9/S9F','fig');saveas(6,'Final_Figures/SupMat/S9/S9F','epsc')
saveas(7,'Final_Figures/SupMat/S9/S9G','fig');saveas(7,'Final_Figures/SupMat/S9/S9G','epsc')
saveas(8,'Final_Figures/SupMat/S9/S9H','fig');saveas(8,'Final_Figures/SupMat/S9/S9H','epsc')
saveas(9,'Final_Figures/SupMat/S9/S9I-Extra','fig');saveas(9,'Final_Figures/SupMat/S9/S9I-Extra','epsc')
saveas(10,'Final_Figures/SupMat/S9/S9J-Extra','fig');saveas(10,'Final_Figures/SupMat/S9/S9J-Extra','epsc')

%% S10

%% S10.1A-D STL1 0.4M NaCl
close all
Make_Panels_For_Figure_3('STL1',0.4);
saveas(31,'Final_Figures/SupMat/S11/S11ABC','fig');saveas(31,'Final_Figures/SupMat/S11/S11ABC','epsc')
saveas(32,'Final_Figures/SupMat/S11/S11D','fig');saveas(32,'Final_Figures/SupMat/S11/S11D','epsc')

% S10.2E-H CTT1 0.2M NaCl
close all
Make_Panels_For_Figure_3('CTT1',0.2);
saveas(31,'Final_Figures/SupMat/S11/S11EFG','fig');saveas(31,'Final_Figures/SupMat/S11/S11EFG','epsc')
saveas(32,'Final_Figures/SupMat/S11/S11H','fig');saveas(32,'Final_Figures/SupMat/S11/S11H','epsc')

% S10.3I-L CTT1 0.4M NaCl
close all
Make_Panels_For_Figure_3('CTT1',0.4);
saveas(31,'Final_Figures/SupMat/S11/S11IJK','fig');saveas(31,'Final_Figures/SupMat/S11/S11IJK','epsc')
saveas(32,'Final_Figures/SupMat/S11/S11L','fig');saveas(32,'Final_Figures/SupMat/S11/S11L','epsc')

%% S12 Metropolis Hastings Convergence
close all
Make_SI_Fig_12
saveas(100,'Final_Figures/SupMat/S12/S12','fig');saveas(100,'Final_Figures/SupMat/S12/S12','epsc')

%% S13/S14 Bias and Uncertainty Plots
close all
Make_Biase_and_Uncertainty_PColors
saveas(12,'Final_Figures/SupMat/S13/S13A','fig');saveas(12,'Final_Figures/SupMat/S13/S13A','epsc')
saveas(17,'Final_Figures/SupMat/S13/S13B','fig');saveas(17,'Final_Figures/SupMat/S13/S13B','epsc')
saveas(13,'Final_Figures/SupMat/S13/S13C','fig');saveas(13,'Final_Figures/SupMat/S13/S13C','epsc')
saveas(18,'Final_Figures/SupMat/S13/S13D','fig');saveas(18,'Final_Figures/SupMat/S13/S13D','epsc')
saveas(14,'Final_Figures/SupMat/S13/S13G','fig');saveas(14,'Final_Figures/SupMat/S13/S13G','epsc')
saveas(19,'Final_Figures/SupMat/S13/S13H','fig');saveas(19,'Final_Figures/SupMat/S13/S13H','epsc')

saveas(2,'Final_Figures/SupMat/S14/S14A','fig');saveas(2,'Final_Figures/SupMat/S14/S14A','epsc')
saveas(7,'Final_Figures/SupMat/S14/S14B','fig');saveas(7,'Final_Figures/SupMat/S14/S14B','epsc')
saveas(3,'Final_Figures/SupMat/S14/S14C','fig');saveas(3,'Final_Figures/SupMat/S14/S14C','epsc')
saveas(8,'Final_Figures/SupMat/S14/S14D','fig');saveas(8,'Final_Figures/SupMat/S14/S14D','epsc')
saveas(4,'Final_Figures/SupMat/S14/S14G','fig');saveas(4,'Final_Figures/SupMat/S14/S14G','epsc')
saveas(9,'Final_Figures/SupMat/S14/S14H','fig');saveas(9,'Final_Figures/SupMat/S14/S14H','epsc')

saveas(5,'Final_Figures/SupMat/S14/S14E','fig');saveas(5,'Final_Figures/4th_Moments/S14E','epsc')
saveas(10,'Final_Figures/SupMat/S14/S14F','fig');saveas(10,'Final_Figures/4th_Moments/S14F','epsc')
saveas(15,'Final_Figures/SupMat/S13/S13E','fig');saveas(15,'Final_Figures/SupMat/S13/S13E','epsc')
saveas(20,'Final_Figures/SupMat/S13/S13F','fig');saveas(20,'Final_Figures/SupMat/S13/S13F','epsc')

saveas(6,'Final_Figures/SupMat/S14/S14I-Extra','fig');saveas(6,'Final_Figures/SupMat/S14/S14I-Extra','epsc')
saveas(11,'Final_Figures/SupMat/S14/S14J-Extra','fig');saveas(11,'Final_Figures/SupMat/S14/S14J-Extra','epsc')
saveas(16,'Final_Figures/SupMat/S13/S13I-Extra','fig');saveas(16,'Final_Figures/SupMat/S13/S13I-Extra','epsc')
saveas(21,'Final_Figures/SupMat/S13/S13J-Extra','fig');saveas(21,'Final_Figures/SupMat/S13/S13J-Extra','epsc')


%% S15 Bias/Uncertainty/TS for CTT1
close all
Make_Panels_For_Figure_4('CTT1')
saveas(41,'Final_Figures/SupMat/S15/S15A','fig');saveas(41,'Final_Figures/SupMat/S15/S15A','epsc')
saveas(42,'Final_Figures/SupMat/S15/S15B','fig');saveas(42,'Final_Figures/SupMat/S15/S15B','epsc')
saveas(43,'Final_Figures/SupMat/S15/S15C','fig');saveas(43,'Final_Figures/SupMat/S15/S15C','epsc')
saveas(44,'Final_Figures/SupMat/S15/S15D','fig');saveas(44,'Final_Figures/SupMat/S15/S15D','epsc')
saveas(45,'Final_Figures/SupMat/S15/S15E','fig');saveas(45,'Final_Figures/SupMat/S15/S15E','epsc')
saveas(46,'Final_Figures/SupMat/S15/S15F','fig');saveas(46,'Final_Figures/SupMat/S15/S15F','epsc')
saveas(47,'Final_Figures/SupMat/S15/S15G','fig');saveas(47,'Final_Figures/SupMat/S15/S15G','epsc')
saveas(48,'Final_Figures/SupMat/S15/S15H','fig');saveas(48,'Final_Figures/SupMat/S15/S15H','epsc')

%% S16 Comparison of FSP and SSA for TS analysis
close all
Compare_FSP_and_SSA_TS
saveas(15,'Final_Figures/SupMat/S16/S16A','fig');saveas(15,'Final_Figures/SupMat/S16/S16A','epsc')
saveas(16,'Final_Figures/SupMat/S16/S16B','fig');saveas(16,'Final_Figures/SupMat/S16/S16B','epsc')
saveas(115,'Final_Figures/SupMat/S16/S16C','fig');saveas(115,'Final_Figures/SupMat/S16/S16C','epsc')
saveas(116,'Final_Figures/SupMat/S16/S16D','fig');saveas(116,'Final_Figures/SupMat/S16/S16D','epsc')
