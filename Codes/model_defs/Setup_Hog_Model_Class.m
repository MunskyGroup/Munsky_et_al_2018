function [ Model ] = Setup_Hog_Model_Class(Model, obj, Parameters)
%Setup_Hog_Model
% This function will take the specified MODEL, the requested
% ANALYSIS and the specified CONDITIONS and PARAMETERS and return the
% resulting Model.PARS for subsequent solution.

%% First to specifiy the Analysis Type:
switch obj.Type
    case {'Moments','Moments_Plus','Moments_Corrected','Moments_4'}
        %     disp('Performing Moment Based Analysis')
        Model.PARS.Moment_Order = obj.Moment_Order; % Order of the moments (=1 or 2).
    case {'Distributions','Distributions_TS'}
    case 'SSA'
    otherwise
        error('Analysis Type Not Recognized');
end

%% Next to specify the Model Type:
% Model
% if strcmp(Model.Type,'N_State_M_Gene')
    Model.PARS.Num_States = Model.NStates; % Number of states for each gene, assuming independent.
    Model.PARS.Num_Genes=Model.MGenes; %  % Number of different genes.
% else
%     error('Model Type Not Recognized');
% end
Model.PARS.Arp8_Gcn5=obj.Hog1p_Type; % Type of Hog1p signal {'WT','ARP8','GCN5','Hot1-5x'}

%% Next to Specify the Conditions:
% Model.PARS.Salt=obj.Salt; % Level of Salt. Only coded for 'WT'.
% Model.PARS.Time_Array=Conditions.Times; % Output times in seconds.

%% Finally to specify how and/or where the parameters are to be defined:
Model.PARS.Pars_Type = Model.Pars_Type; % Type of parameters.
if strcmp(Model.PARS.Pars_Type,'N_State_General')||strcmp(Model.PARS.Pars_Type,'N_State_General_OLD_Format');
    Load_File_Names = Parameters.Load_File_Names; % Names of parameter files.
    for i=1:Model.PARS.Num_Genes
        load(Load_File_Names{i});
        Model.PARS.Gene(i).PARS = BEST_PARS_ALL; % Set parameters of the first gene.
    end
    if strcmp(Model.PARS.Pars_Type,'N_State_General_OLD_Format');
        Model.PARS = Reformat_PARS_Definition(Model.PARS);
        Model.PARS.Pars_Type = 'N_State_General';
    end
    Model.PARS = Reformat_PARS_Definition_to_Vector(Model.PARS);
    BEST_PARS_VECTOR = Model.PARS.Vector;
    save('Parameters_Joint/PARS_Vector.mat','BEST_PARS_VECTOR');
elseif strcmp(Model.PARS.Pars_Type,'Vector')
    load(Parameters.Load_File_Name); % Load parameter files.
    Model.PARS.Vector = BEST_PARS;
    Model.PARS = Reformat_PARS_Definition_Vector_Class(Model.PARS);
    Model.PARS.Pars_Type = 'N_State_General';
elseif strcmp(Model.PARS.Pars_Type,'Vector_Given')
    Model.PARS.Vector = Parameters.Vector;
    Model.PARS = Reformat_PARS_Definition_Vector_Class(Model.PARS);
    Model.PARS.Pars_Type = 'N_State_General';
elseif strcmp(Model.PARS.Pars_Type,'Connected')
    Model.PARS.Vector = Parameters.Vector;
    Model.PARS.Model_Connection = Model.Connected;
    Model.PARS = PARS_Definition_Connected_Models(Model.PARS);
    Model.PARS.Pars_Type = 'Connected';
else
    error('Parameter Type Not Recognized');
end

if strcmp(obj.Spatial,'Fixed')
    Model.PARS.Degrade(:,2)=Model.PARS.Degrade(:,1);  % only one degradation rate for the cytoplasmic RNA.
end
