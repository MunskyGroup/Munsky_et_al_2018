function Reset_GUI(GUI)
GUI.MOD.Mod_Ovrd_Struct = GUI.MOD.Model_Obj;
GUI.MOD.Salt = [0.2 0.4];

GUI.MOD.Model_Override = 1;

GUI.MOD.Mod_Ovrd_Struct.PARS.States(1,1,2) = GUI.k12EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.States(1,2,3) = GUI.k23EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.States(1,3,4) = GUI.k34EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.States(1,3,2) = GUI.k32EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.States(1,4,3) = GUI.k43EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.States(1,2,1) = GUI.k21aEditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.States_Time(1,2,1) = GUI.k21bEditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.Prod(1) = GUI.kr1EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.Prod(2) = GUI.kr2EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.Prod(3) = GUI.kr3EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.Prod(4) = GUI.kr4EditField.Value;
GUI.MOD.Mod_Ovrd_Struct.PARS.t_offset = GUI.timeoffsetEditField.Value;
switch GUI.MOD.Spatial
    case 'Spatial'
        GUI.MOD.Mod_Ovrd_Struct.PARS.Degrade(2) = GUI.degradationEditField.Value;
        GUI.MOD.Mod_Ovrd_Struct.PARS.Degrade(1) = GUI.nucleardegradationEditField.Value;
        GUI.MOD.Mod_Ovrd_Struct.PARS.Transport = GUI.transportEditField.Value;
    case 'NonSpatial'
        GUI.MOD.Mod_Ovrd_Struct.PARS.Degrade(1) = GUI.degradationEditField.Value;
end

GUI.Moments = [];
GUI.Distributions = [];

