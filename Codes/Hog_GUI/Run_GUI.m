close all force
clear all
GUI = Hog_Gui;

%%
path(path,'../../Munsky_2017_Final_Hog_Codes')
path(path,'../Plotting_Functions');
path(path,'../Hog_Data');
MOD = Hog_Model;

MOD.Gene = GUI.GeneTypeButtonGroup.SelectedObject.Text;
disp(['Gene = ',MOD.Gene,'.']);

switch GUI.AnalysisTypeButtonGroup.SelectedObject.Text
    case 'FSP'
        disp('Parameters found with FSP');
        MOD.Fit_Type='Distributions';
    case 'Means'
        disp('Parameters found with Means');
        MOD.Fit_Type='Means';
    case '2nd Moments'
        disp('Parameters found with Means and Variances');
        MOD.Fit_Type='Moments_Corrected';
    case 'Extended Moments'
        disp('Parameters found with Extended Moments');
        MOD.Fit_Type='Moments_4';
end

MOD.Spatial = GUI.SpatialTypeButtonGroup.SelectedObject.Text;
disp(['Using ',MOD.Spatial,' parameters.']);

DIR = ['../Fit_Results/',MOD.Gene,'_0p2and0p4/',MOD.Spatial,'/',MOD.Fit_Type,'/',MOD.Model];
MOD.Parameter_File = [DIR,'/Best_Pars_AC.mat'];

