function Moments_Out = Compute_Moment_Trajectories(Moment_Properties)
% This function will take a structure input that defines a linear model,
% any time varying components, and the desired output charateristics. 
% Properties of the input should include:
% Moment_Properties.Order = {1,[2]} -- Moment order.
% Moment_Properties.Output_Times = [linspace(0,1,10)] -- times to compute outputs.
% Moment_Properties = propensity function structure.
%   Moment_Properties.Stoichiometry = Stoichiometry Matrix
%   Moment_Properties.W0 = zeroth order propensities.
%   Moment_Properties.W1 = first order propensities.
% Moment_Properties.Signal = time varying signal structure
%   Moment_Properties.Signal.Type = {'None','Interpolate',['Function']}  -- Type of time varying signal.
%   Moment_Properties.Signal.Function = @(t)Define_Hog_Signal(t,Reactions) -- Function handle if time varying function.
%   Moment_Properties.Signal.Time -- Vector of times for empiricle signal if "interpolation"
%   Moment_Properties.Signal.Value -- Vector of Values for empiricle signal if "interpolation"
%   Moment_Properties.Signal.W0t = []; zeroth order time varying propensities.
%   Moment_Properties.Signal.W1t = []; first order time varying propensities.
% Moment_Properties.Moments_0 -- Initial moments at first time point.    
% Moment_Properties.N_species -- Number of species (number of rows in Stoichometry).    

Jacobian_Fun = @(t,x)Get_Moments_Jacobian(t,[],Moment_Properties.Stoichiometry,...
    Moment_Properties.W0,Moment_Properties.W1,Moment_Properties.Moment_Order,[],Moment_Properties.Signal);

Moment_Properties.Signal.Jac = Get_Moment_Jac(Moment_Properties.Stoichiometry,...
    Moment_Properties.Signal.W0t,Moment_Properties.Signal.W1t,Moment_Properties.Moment_Order);
Moment_Properties.Jac = Get_Moment_Jac(Moment_Properties.Stoichiometry,...
    Moment_Properties.W0,Moment_Properties.W1,Moment_Properties.Moment_Order);

W0_SAT = Moment_Properties.W0 + Moment_Properties.Signal.W0t*1e10;   % Saturation of W0
W1_SAT = Moment_Properties.W1 + Moment_Properties.Signal.W1t*1e10;   % Saturation of W1
W0_SAT = sparse(max(0,W0_SAT));  %Propensities must always be positive for any State.
W1_SAT = sparse(max(0,W1_SAT));  %Propensities must always be positive for any State.
Moment_Properties.Jac_SAT = Get_Moment_Jac(Moment_Properties.Stoichiometry,W0_SAT,W1_SAT,Moment_Properties.Moment_Order);

JPattern = sparse(((Moment_Properties.Signal.Jac+Moment_Properties.Jac)~=0));

ode_options = odeset('MaxStep',10,'Jacobian',Jacobian_Fun,'JPattern',JPattern);

fun = @(t,Moments)Get_Moments_Derivative_JacBased(t,Moments,...
    Moment_Properties.Stoichiometry,...
    Moment_Properties.W0,...
    Moment_Properties.W1,...
    Moment_Properties.N_species,...
    Moment_Properties.Moment_Order,...
    Moment_Properties.Signal,...
    Moment_Properties.Jac,...
    Moment_Properties.Signal.Jac,...
    Moment_Properties.Jac_SAT);

[Moments_Out.Times,Moments_Out.Trajectories] = ode23s(fun,...
    Moment_Properties.Output_Times,...
    Moment_Properties.Moments_0,...
    ode_options);  %solve CME at given initiation times.