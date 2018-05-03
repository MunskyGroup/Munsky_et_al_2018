function Rest_GUI_Parameters(GUI)
%% Update Parameter Tables
path(path,'Plotting_Functions');
path(path,'Hog_Data');
GUI.MOD.Model_Override=0;
GUI.MOD.Gene = GUI.GeneTypeButtonGroup.SelectedObject.Text;
disp(['Gene = ',GUI.MOD.Gene,'.']);

switch GUI.AnalysisTypeButtonGroup.SelectedObject.Text
    case 'FSP'
        disp('Parameters found with FSP');
        GUI.MOD.Fit_Type='Distributions';
    case 'Means'
        disp('Parameters found with Means');
        GUI.MOD.Fit_Type='Means';
    case '2nd Moments'
        disp('Parameters found with Means and Variances');
        GUI.MOD.Fit_Type='Moments_Corrected';
    case 'Extended Moments'
        disp('Parameters found with Extended Moments');
        GUI.MOD.Fit_Type='Moments_4';
end

GUI.MOD.Spatial = GUI.SpatialTypeButtonGroup.SelectedObject.Text;
disp(['Using ',GUI.MOD.Spatial,' parameters.']);

DIR = ['Fit_Results/',GUI.MOD.Gene,'_0p2and0p4/',GUI.MOD.Spatial,'/',GUI.MOD.Fit_Type,'/',GUI.MOD.Model];
GUI.MOD.Parameter_File = [DIR,'/Best_Pars_AC.mat'];

GUI.k12EditField.Value = GUI.MOD.Model_Obj.PARS.States(1,1,2);
GUI.k23EditField.Value = GUI.MOD.Model_Obj.PARS.States(1,2,3);
GUI.k34EditField.Value = GUI.MOD.Model_Obj.PARS.States(1,3,4);
GUI.k32EditField.Value = GUI.MOD.Model_Obj.PARS.States(1,3,2);
GUI.k43EditField.Value = GUI.MOD.Model_Obj.PARS.States(1,4,3);
GUI.k21aEditField.Value = GUI.MOD.Model_Obj.PARS.States(1,2,1);
GUI.k21bEditField.Value = GUI.MOD.Model_Obj.PARS.States_Time(1,2,1);
GUI.kr1EditField.Value = GUI.MOD.Model_Obj.PARS.Prod(1);
GUI.kr2EditField.Value = GUI.MOD.Model_Obj.PARS.Prod(2);
GUI.kr3EditField.Value = GUI.MOD.Model_Obj.PARS.Prod(3);
GUI.kr4EditField.Value = GUI.MOD.Model_Obj.PARS.Prod(4);
GUI.timeoffsetEditField.Value = GUI.MOD.Model_Obj.PARS.t_offset;
switch GUI.MOD.Spatial
    case 'Spatial'
        GUI.degradationEditField.Value = GUI.MOD.Model_Obj.PARS.Degrade(2);
        GUI.nucleardegradationEditField.Value = GUI.MOD.Model_Obj.PARS.Degrade(1);
        GUI.transportEditField.Value = GUI.MOD.Model_Obj.PARS.Transport;
    case 'NonSpatial'
        GUI.degradationEditField.Value = GUI.MOD.Model_Obj.PARS.Degrade(1);
        GUI.nucleardegradationEditField.Value = 0;
end
GUI.mRNAelongationEditField.Value = 63;
GUI.Moments=[];
GUI.Distributions=[];
