function [Dynamics,N_Space] = hog_model_FSP_class(obj)
% This function will use the parameter set (PARS) anf the number of states,
% and the mutatn type (Hog1p_Type), salt level, and number of genes, and
% then run the program to compute the first and second moments as a
% function of time.
%% N RNA species, independent four states (4^N State Total)
PARS = obj.Model_Obj.PARS;
if obj.MGenes~=1&&(strcmp(obj.Spatial,'Spatial')||strcmp(obj.Spatial,'Fixed'))
    error('Not set up yet to handle multiple genes and spatial analyses for FSP analysis')
end

if strcmp(PARS.Pars_Type,'Independent_OLD_Format');
    error('FSP not ready for this pars type yet')
    %     PARS.Pars_independent = []; % Total set of Parameters, assuming that all mRNA
    %     % products have independent states.
    %     for i=1:obj.MGenes
    %         PARS.Pars_independent = [PARS.Pars_independent;PARS.Gene(i).PARS];
    %     end
    %     Model.G_States = obj.Num_States^obj.MGenes; %total number of states.
    %     [Reactions,Model.Initial_Condition] = Define_Pars_NGenes(PARS,obj.MGenes,obj.Spatial);
    %     % Call function to Define the Linear ODE model for the Moment closure,
    %     % including the Reactions of the process and the initial condition for the
    %     % moments.  The Output 'Reactions' is a structure with the substructures:
    %     % Reactions.Hog_Parameters - Paramters of the empiricle Hog1p model,
    %     % Reactions.Pars_State_Change - Parameters for all of the gene transition reactions
    %     % Reactions.Pars_Prods - Parameter for the Production events,
    %     % Reactions.Pars_Degrade_Nuc - Nuclear mRNA degradation parameters.
    %     % Reactions.Pars_Degrade_CYt - Cytoplasmic mRNA degradation parameters.
    %     % Reactions.Pars_Transport - mRNA transport rate
    
elseif strcmp(PARS.Pars_Type,'N_State_General');
    Model.G_States = obj.Num_States^obj.MGenes; %total number of states.
    [Reactions] = Define_Pars_NGenes(obj,PARS);
elseif strcmp(PARS.Pars_Type,'Connected');
    %     error('FSP not ready for this pars type yet')
    Model.G_States = obj.Num_States; %total number of states.
    [Reactions] = Define_Pars_Connected(obj,PARS);
end

Model.Signal.Type = 'Function'; % Type of input signal {'Function','
if strcmp(obj.Hog_Signal,'Default')
    Model.Signal.Function = @(t)Define_Hog_Signal(t,Reactions.Hog_Parameters);
else 
    Model.Signal.Function = Analysis.Hog_Signal_Function;
end
% Call subroutine to define the Hog1p signal using the parameters above.

% Formulate structure to define the model for the integration.
Model.Pars_State_Change = Reactions.Pars_State_Change;
Model.Pars_Prods = Reactions.Pars_Prods;
Model.Pars_Transport = Reactions.Pars_Transport;
Model.Pars_Degrade_Nuc = Reactions.Pars_Degrade_Nuc;
Model.Pars_Degrade_Cyt = Reactions.Pars_Degrade_Cyt;

if strcmp(obj.Spatial,'Fixed')
    Model.Pars_Degrade_Cyt = Model.Pars_Degrade_Nuc;
    Reactions.Pars_Degrade_Cyt = Reactions.Pars_Degrade_Nuc;
end

PARS.HISTOGRAMS = @Get_HIST_General_Spatial;
% PARS.CME = @cme_General_Spatial;
PARS.obj.Num_States = obj.Num_States;

if strcmp(obj.Spatial,'Spatial')||strcmp(obj.Spatial,'Fixed')&&obj.MGenes==1
    PARS.N_Space = [obj.Num_States 15 200];
elseif strcmp(obj.Spatial,'NonSpatial')&&obj.MGenes==1
    PARS.N_Space = [obj.Num_States 300 0];
elseif strcmp(obj.Spatial,'NonSpatial')&&obj.MGenes==2&&strcmp(Analysis.Gene,'STL1')
    PARS.N_Space = [obj.Num_States 300 0];    
elseif strcmp(obj.Spatial,'NonSpatial')&&obj.MGenes==2&&strcmp(Analysis.Gene,'CTT1')
    PARS.N_Space = [obj.Num_States 0 300];    
else
    error('The FSP is not ready for this many indepedent species');    
end
PARS.N_states=PARS.N_Space(1)*...
    (PARS.N_Space(2)+1)*...
    (PARS.N_Space(3)+1);  % Maximum number of states.

if ~isfield(PARS,'JPattern')||isempty(PARS.JPattern)
    PARS.Spatial=obj.Spatial;
    PARS = cme_General_Spatial_NState(PARS);  %% This function will define and store the infintesimal generator of the master equation.
end
PARS.Time_Array = obj.tt*60;
Dynamics = PARS.HISTOGRAMS(PARS,Model,Reactions.Hog_Parameters);
% Call sub routine to solve the model for the moments at the specified time
% points.
N_Space=PARS.N_Space;

function [Reactions] = Define_Pars_NGenes(obj,PARS)
% This function defines all of the reactions and parameters for the Hog1p
% activated transcription of N different independent genes.
% Required inputs are:
% The parameters in the structure 'PARS',
% 'Hog1p_Type = {'WT','ARP8','GCN5','Hot1-5x'} -- The type of syste that
% generates the Hog1p signal.
% 'obj.Num_States' -- the number of gene states per mRNA product.
% 'SALT' the amount of NaCl (only for 'WT')
% 'N_Prods' -- the number of mRNA products
%
% The Output 'Reactions' is a structure with the substructures:
% Reactions.Hog_Parameters - Paramters of the empiricle Hog1p model,
% Reactions.Pars_State_Change - Parameters for all of the gene transition reactions
% Reactions.Pars_Prods - Parameter for the Production events,
% Reactions.Pars_Degrade_Nuc - Nuclear mRNA degradation parameters.
% Reactions.Pars_Degrade_CYt - Cytoplasmic mRNA degradation parameters.
% Reactions.Pars_Transport - mRNA transport rateReactions=[];

if strcmp(obj.Hog1p_Type,'WT') % Wild Type    
    Reactions.Hog_Parameters.eta =5.890086837674720e+00;
    Reactions.Hog_Parameters.r1 = 6.094509227679889e-03;
    alpha = 1.679176848801883e-03;
    Reactions.Hog_Parameters.M = 4.651137270724188e+01;
    Reactions.Hog_Parameters.del = -4.096887973322170e-02;
    Reactions.Hog_Parameters.t0 = 5.509466971325969e+01;
    Reactions.Hog_Parameters.r2 = alpha/(obj.Salt-Reactions.Hog_Parameters.del);
elseif strcmp(obj.Hog1p_Type,'ARP8') % Arp8 Mutant
    error('Hog1p Signal for ARP8 not up to date')
    Reactions.Hog_Parameters.eta =1.56653632839433;
    Reactions.Hog_Parameters.r1 = 2.23892760575769e-09;
    Reactions.Hog_Parameters.M = 2458185.38152383;
    Reactions.Hog_Parameters.r2 = 0.00257402402551422;
elseif strcmp(obj.Hog1p_Type,'GCN5') % Gcn5 Mutant
    error('Hog1p Signal for GCN5 not up to date')
    Reactions.Hog_Parameters.eta =0.79163;
    Reactions.Hog_Parameters.r1 = 4.9453e-16;
    Reactions.Hog_Parameters.M = 2.8321e+12;
    Reactions.Hog_Parameters.r2 = 0.0032722;
elseif strcmp(obj.Hog1p_Type,'Hot1-5x') % Hot1-5x Mutant
    error('Hog1p Signal for Hot1-5x not up to date')
    Reactions.Hog_Parameters.eta =0.87595;
    Reactions.Hog_Parameters.r1 = 2.8981e-14;
    Reactions.Hog_Parameters.M = 3.1105e+11;
    Reactions.Hog_Parameters.r2 = 0.0034319;
end
if strcmp(PARS.Pars_Type,'Independent_OLD_Format');
    %     switch PARS.obj.Num_States
    %         case 2
    %         case 3
    %         case 4
    %             Reactions = Assign_Reactions_4State(Reactions,PARS,N_Prods);
    %             Reactions.Hog_Parameters
    %             PARS.Gene(1).PARS(25)
    %             PARS.Pars_independent
    %         case 5
    %     end
else
    Reactions.Hog_Parameters.t_offset=PARS.t_offset;
    Reactions = assign_reactions_nstate_class(Reactions,obj,PARS);
end

%% function [Reactions,Initial_Condition] = Define_Pars_Connected(PARS,N_Prods,Spatial)
function [Reactions] = Define_Pars_Connected(PARS,N_Prods,Spatial)
% This function defines all of the reactions and parameters for the Hog1p
% activated transcription of N (=2) different correlated genes.
% Required inputs are:
% The parameters in the structure 'PARS',
% 'Hog1p_Type = {'WT','ARP8','GCN5','Hot1-5x'} -- The type of syste that
% generates the Hog1p signal.
% 'obj.Num_States' -- the number of gene states per mRNA product.
% 'SALT' the amount of NaCl (only for 'WT')
% 'N_Prods' -- the number of mRNA products
%
% The Output 'Reactions' is a structure with the substructures:
% Reactions.Hog_Parameters - Paramters of the empiricle Hog1p model,
% Reactions.Pars_State_Change - Parameters for all of the gene transition reactions
% Reactions.Pars_Prods - Parameter for the Production events,
% Reactions.Pars_Degrade_Nuc - Nuclear mRNA degradation parameters.
% Reactions.Pars_Degrade_CYt - Cytoplasmic mRNA degradation parameters.
% Reactions.Pars_Transport - mRNA transport rateReactions=[];

if strcmp(obj.Hog1p_Type,'WT') % Wild Type
    Reactions.Hog_Parameters.eta =5.890086837674720e+00;
    Reactions.Hog_Parameters.r1 = 6.094509227679889e-03;
    alpha = 1.679176848801883e-03;
    Reactions.Hog_Parameters.M = 4.651137270724188e+01;
    Reactions.Hog_Parameters.del = -4.096887973322170e-02;
    Reactions.Hog_Parameters.t0 = 5.509466971325969e+01;
    Reactions.Hog_Parameters.r2 = alpha/(obj.Salt-Reactions.Hog_Parameters.del);
elseif strcmp(obj.Hog1p_Type,'ARP8') % Arp8 Mutant
    Reactions.Hog_Parameters.eta =1.56653632839433;
    Reactions.Hog_Parameters.r1 = 2.23892760575769e-09;
    Reactions.Hog_Parameters.M = 2458185.38152383;
    Reactions.Hog_Parameters.r2 = 0.00257402402551422;
elseif strcmp(obj.Hog1p_Type,'GCN5') % Gcn5 Mutant
    Reactions.Hog_Parameters.eta =0.79163;
    Reactions.Hog_Parameters.r1 = 4.9453e-16;
    Reactions.Hog_Parameters.M = 2.8321e+12;
    Reactions.Hog_Parameters.r2 = 0.0032722;
elseif strcmp(obj.Hog1p_Type,'Hot1-5x') % Hot1-5x Mutant
    Reactions.Hog_Parameters.eta =0.87595;
    Reactions.Hog_Parameters.r1 = 2.8981e-14;
    Reactions.Hog_Parameters.M = 3.1105e+11;
    Reactions.Hog_Parameters.r2 = 0.0034319;
end

Reactions.Hog_Parameters.t_offset=PARS.t_offset;
Reactions = Assign_Reactions_Connected(Reactions,PARS,N_Prods);
% 
% if strcmp(Spatial,'yes')
%     Initial_Condition.Means=zeros(PARS.obj.Num_States+N_Prods*2,1);
%     Initial_Condition.Means(1)=1; % Start in State 1 with no RNA.
%     Initial_Condition.SIG2 = zeros(PARS.obj.Num_States+N_Prods*2,PARS.obj.Num_States+N_Prods*2); % Delta distribution.
% else
%     Initial_Condition.Means=zeros(PARS.obj.Num_States+N_Prods,1);
%     Initial_Condition.Means(1)=1; % Start in State 1 with no RNA.
%     Initial_Condition.SIG2 = zeros(PARS.obj.Num_States+N_Prods,PARS.obj.Num_States+N_Prods); % Delta distribution.
% end