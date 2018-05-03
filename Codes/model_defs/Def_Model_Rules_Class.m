function Model = Def_Model_Rules_Class(obj)%% Apply_Model_Rules
%This function define the rules for specifying the model parameters.  This
%is where we need to make changes when we add new models.  For example, if
%we want to remove certain hog1-dependent reactions.  Or if we want to
%define the file for parameter sets fro a different analysis type. 

Model.Pos_Pars = [];  
%These are the parameters that are known to be positive.  For best
%efficiency, we should define these according to the sp[ecific model
%specification.

switch obj.Model
    case '4State_k21only_noloop'
        Model.MGenes=1; % Number of different genes.
        Model.Connected = [];
        Model.Free_Pars = [1 4 5 8 9 12,12+4,25:28,33,39];
        Model.Pos_Pars = [1 5 8 9 12,25:28,33,39];
        if strcmp(obj.Spatial,'Spatial')
            Model.Free_Pars = [Model.Free_Pars,35,37];  % Add deg in Cyt and Transport.
            Model.Pos_Pars = [Model.Pos_Pars,35,37];  % Add deg in Cyt and Transport.
        elseif strcmp(obj.Spatial,'Fixed')
            Model.Free_Pars = [Model.Free_Pars,37];  % Add Transport.
            Model.Pos_Pars = [Model.Pos_Pars,37];    % Add Transport.
        end
        Model.NStates =4; % Number of states for each gene, assuming independent.
end

if isempty(obj.Parameter_File)
    if obj.MGenes==1
        Model.Base_Par_File_Name = ['Fit_Results/',...
            obj.Gene,'_0p2and0p4/',...
            obj.Spatial,'/',...
            obj.Type,'/',...
            obj.Model,'/Best_Pars',obj.Reps,'.mat'];
    else
        Model.Base_Par_File_Name = ['Fit_Results/',...
            obj.Gene,'/',...
            obj.Spatial,'/',...
            obj.Type,'/',...
            obj.Model,'/Best_Pars',obj.Reps,'.mat']; 
    end
else
    Model.Base_Par_File_Name = obj.Parameter_File;
end

switch obj.Model
    case {'4State_k21only_noloop'}
        Model.Pars_Type = 'Vector';  % Here the model parameters are stored in a vector in an m-file.
            Parameters.Load_File_Name= Model.Base_Par_File_Name;
        Model = Setup_Hog_Model_Class(Model, obj, Parameters);
        Model.Pars_Type = 'Vector_Given';  % Here the model parameters are stored in a vector in an m-file.
end