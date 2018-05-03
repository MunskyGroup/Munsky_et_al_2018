function [Model_Properties,Moments_Out]=define_stoich_class(Model,obj)
% Input structure, 'Model', includes fields:
% Model.G_States = Number of phenotype gene states.
% obj.MGenes = Number of product species
% Model.Spatial = {'no',['yes']} Is this a spatial model.
% Model.Pars_State_Change.Const -- Matrix of parameters for state changes.
% Model.Pars_State_Change.Time -- Matrix of time varying parameters for state changes.
% Model.Pars_Prods.Const -- matrix of parameters for production rates (each
%       row is for a different product species.  Each column is for a differen
%       phenotype state.
% Model.Pars_Prods.Time -- matrix of parameters for time varying production rates (each
%       row is for a differen product species.  Each column is for a differen
%       phenotype state.
% Model.Pars_Transport.Const -- vector of parameters for transport rates (each
%       column is for a different product species.
% Model.Pars_Transport.Time -- vector of parameters for time varying transport rates (each
%       column is for a different product species.
% Model.Pars_Degrade_Nuc.Const -- vector of parameters for degradation rates (each
%       column is for a different product species.
% Model.Pars_Degrade_Nuc.Time -- vector of parameters for time varying degradation rates (each
%       column is for a different product species.
% Model.Pars_Degrade_Cyt.Const -- vector of parameters for degradation rates (each
%       column is for a different product species.
% Model.Pars_Degrade_Cyt.Time -- vector of parameters for time varying degradation rates (each
%       column is for a different product species.
% Model.Initial_Condition.Means -- Initial mean levels.
% Model.Initial_Condition.SIG2 -- Initial covariance matrix
% Model.Time_Array -- Time at which to solve for model properties

switch obj.Spatial
    case {'Spatial','Fixed'}
        N_Species = Model.G_States+2*obj.MGenes;  % The total number of species is the number of states plus nuclear and cytoplasmic products.
        N_reactions = sum(sum((abs(Model.Pars_State_Change.Const)+abs(Model.Pars_State_Change.Time))~=0))+... % State Changes.
            obj.MGenes*4; % Production, transport, and Degradation (nuc and cyt).
        % Total number of non-zero reactions.
    case 'NonSpatial'
        N_Species = Model.G_States+obj.MGenes;  % The total number of species is the number of states plus mRNA products.
        N_reactions = sum(sum((abs(Model.Pars_State_Change.Const)+abs(Model.Pars_State_Change.Time))~=0))+... % State Changes.
            obj.MGenes*2; % Production and Degradation.
        % Total number of non-zero reactions.
end

%% Initialize Stoichiometry and Propensity functions.

Model_Properties.Stoichiometry = spalloc(N_Species,N_reactions,N_reactions*2);  % Stoichiometry matrix
Model_Properties.W0 = zeros(N_reactions,1); % Zeroth order reactions
Model_Properties.W1 = zeros(N_reactions,N_Species); % First order reactions
Model_Properties.Signal = Model.Signal;
Model_Properties.Signal.W0t = zeros(N_reactions,1); % Zeroth order reactions (Time Varying)
Model_Properties.Signal.W1t = zeros(N_reactions,N_Species);% First order reactions (Time Varying)

%% Define Stoichiometry and Rates for Phenotype Changes.
k=0;
for i=1:Model.G_States
    for j=i+1:Model.G_States
        if (Model.Pars_State_Change.Const(i,j)~=0)||(Model.Pars_State_Change.Time(i,j)~=0)
            k=k+1;
            Model_Properties.Stoichiometry([i,j],k) = [-1 1]';  % Si -> Sj
            Model_Properties.W1(k,i) = Model.Pars_State_Change.Const(i,j);
            Model_Properties.Signal.W1t(k,i) = Model.Pars_State_Change.Time(i,j);
        end
        if (Model.Pars_State_Change.Const(j,i)~=0)||(Model.Pars_State_Change.Time(j,i)~=0)
            k=k+1;
            Model_Properties.Stoichiometry([i,j],k) = [1 -1]';  % Sj -> Si
            Model_Properties.W1(k,j) = Model.Pars_State_Change.Const(j,i);
            Model_Properties.Signal.W1t(k,j) = Model.Pars_State_Change.Time(j,i);
        end
    end
end

%% Define Stoichiometry and Rates for Production.
Model_Properties.transcription_inds=[];
for i=1:obj.MGenes
    k=k+1;
    switch obj.Spatial
        case {'Spatial','Fixed'}
            Model_Properties.Stoichiometry(Model.G_States+(i-1)*2+1,k)=1;  % Productions of RNA
        case 'NonSpatial'
            Model_Properties.Stoichiometry(Model.G_States+i,k)=1;  % Productions of RNA
    end
    %     size(Model_Properties.W1)
    %     size(Model.Pars_Prods.Const)
    %     Model.G_States
    Model_Properties.W1(k,1:Model.G_States) = Model.Pars_Prods.Const(i,:);
    Model_Properties.Signal.W1t(k,1:Model.G_States) = Model.Pars_Prods.Time(i,:);
    Model_Properties.transcription_inds = [Model_Properties.transcription_inds,k];
end

%% Define Stoichiometry and Rates for Degradation of Nuclear RNA
for i=1:obj.MGenes
    k=k+1;
    switch obj.Spatial
        case {'Spatial','Fixed'}
            Model_Properties.Stoichiometry(Model.G_States+(i-1)*2+1,k) =-1;  % Degradation of Nuclear RNA
            Model_Properties.W1(k,Model.G_States+(i-1)*2+1) = Model.Pars_Degrade_Nuc.Const(i);
            Model_Properties.Signal.W1t(k,Model.G_States+(i-1)*2+1) = Model.Pars_Degrade_Nuc.Time(i);
        case 'NonSpatial'
            Model_Properties.Stoichiometry(Model.G_States+i,k) =-1;  % Degradation of Nuclear RNA
            Model_Properties.W1(k,Model.G_States+i) = Model.Pars_Degrade_Nuc.Const(i);
            Model_Properties.Signal.W1t(k,Model.G_States+i) = Model.Pars_Degrade_Nuc.Time(i);
    end
end

if strcmp(obj.Spatial,'Spatial')||strcmp(obj.Spatial,'Fixed')
    %% Define Stoichiometry and Rates for Transport of RNA
    for i=1:obj.MGenes
        k=k+1;
        Model_Properties.Stoichiometry(Model.G_States+(i-1)*2+[1,2],k) =[-1 1]';  % Transport of RNA
        Model_Properties.W1(k,Model.G_States+(i-1)*2+1) = Model.Pars_Transport.Const(i);
        Model_Properties.Signal.W1t(k,Model.G_States+(i-1)*2+1) = Model.Pars_Transport.Time(i);
    end
end

%% Define Stoichiometry and Rates for Degradation of Cytoplasmic RNA
if strcmp(obj.Spatial,'Spatial')
    for i=1:obj.MGenes
        k=k+1;
        Model_Properties.Stoichiometry(Model.G_States+(i-1)*2+2,k) =-1;  % Degradation of Cytoplasmic RNA
        Model_Properties.W1(k,Model.G_States+(i-1)*2+2) = Model.Pars_Degrade_Cyt.Const(i);
        Model_Properties.Signal.W1t(k,Model.G_States+(i-1)*2+2) = Model.Pars_Degrade_Cyt.Time(i);
    end
elseif strcmp(obj.Spatial,'Fixed')
    for i=1:obj.MGenes
        k=k+1;
        Model_Properties.Stoichiometry(Model.G_States+(i-1)*2+2,k) =-1;  % Degradation of Cytoplasmic RNA -- SAME as NUCLEAR.
        Model_Properties.W1(k,Model.G_States+(i-1)*2+2) = Model.Pars_Degrade_Nuc.Const(i);
        Model_Properties.Signal.W1t(k,Model.G_States+(i-1)*2+2) = Model.Pars_Degrade_Nuc.Time(i);
    end
end
%% Define the initial conditions.
try
    Model_Properties.Moments_0 = Model.Initial_Condition.Means;
    if obj.Moment_Order >= 2
        SIG2 = Model.Initial_Condition.SIG2;
        k=N_Species;
        for i=1:N_Species
            for j=i:N_Species
                k=k+1;
                Model_Properties.Moments_0(k,1)=SIG2(i,j);
            end
        end
    end
    if obj.Moment_Order >= 3
        SIG3 = Model.Initial_Condition.SIG3;
        for i1=1:N_Species
            for i2=i1:N_Species
                for i3=i2:N_Species
                    k=k+1;
                    Model_Properties.Moments_0(k,1)=SIG3(i1,i2,i3);
                end
            end
        end
    end
    if obj.Moment_Order>=4
        SIG4 = Model.Initial_Condition.SIG4;
        for i1=1:N_Species
            for i2=i1:N_Species
                for i3=i2:N_Species
                    for i4=i3:N_Species
                        k=k+1;
                        Model_Properties.Moments_0(k,1)=SIG4(i1,i2,i3,i4);
                    end
                end
            end
        end
    end
catch
    disp('Error in Moments IC Specification')
end
%%
Model_Properties.Output_Times = Model.T_array;
Model_Properties.N_species = N_Species;

Model_Properties.Stoichiometry = sparse(Model_Properties.Stoichiometry);
Model_Properties.W0 = sparse(Model_Properties.W0);
Model_Properties.W1 = sparse(Model_Properties.W1);
Model_Properties.Signal.W0t = sparse(Model_Properties.Signal.W0t);
Model_Properties.Signal.W1t = sparse(Model_Properties.Signal.W1t);
Model_Properties.Moment_Order = obj.Moment_Order;

try
    Moments_Out = Compute_Moment_Trajectories(Model_Properties);
catch
    Moments_Out=[];
    disp('Error in moments analysis')
end
% Moments = Moments_Out.Trajectories;