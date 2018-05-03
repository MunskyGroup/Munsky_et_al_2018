function Moments_Jacobian = Get_Moments_Jacobian(t,~,S,W0,W1,Moment_Order,~,Signal)
% This function computes the derivative with respect to time for the
% specified stoichiometry and propensity functions.  At present it will
% only produce the derivatives for the first and second moments. It is
% assumed that the parameter can be constant or varying as a function of
% time. This function can either be an inline function of time or a set of
% data from which the signal can be interpolated.
%
% t = % current time
% Moments = First (and second moments) in array form.
% S =   % Stoichiometry
% W0 =  % Constant (wrt X) Reaction Rates, costant wrt time.
% W1 =  % Linear (wrt X) Reaction Rates, costant wrt time.
% W0t =  % Constant (wrt X) Reaction Rates, linear wrt time.
% W1t =  % Linear (wrt X) Reaction Rates, linear wrt time.
% Signal.Time = time points for the input signal from whoich to
% interpolate.
% Signal.Value = Value points for the input signal from which to
% interpolate.

%%  Find W0 and W1 at current time
if nargin>=7
    if strcmp(Signal.Type,'Interpolate')==1
    it =1;
    while t>Signal.Time(it)&&it<length(Signal.Time)
        it=it+1;
    end
    if it==1
        Input_Signal = Signal.Value(1);
    elseif it<=length(Signal.Time)
        Input_Signal = Signal.Value(it-1)+...
            (Signal.Value(it)-Signal.Value(it-1))*...
            (t-Signal.Time(it-1))/(Signal.Time(it)-Signal.Time(it-1));
    else
        Input_Signal = Signal.Value(length(Mean_ERK));
    end

    elseif strcmp(Signal.Type,'Function')==1
        Input_Signal = Signal.Function(t);
    end
    W0 = W0 + Signal.W0t*Input_Signal;   % Time varying values of W0
    W1 = W1 + Signal.W1t*Input_Signal;   % Time varying values of W1
end
%% Enforce positivity on all rates.
W0 = max(0,W0);  %Propensities must always be positive for any State.
W1 = max(0,W1);  %Propensities must always be positive for any State.

Moments_Jacobian = Get_Moment_Jac(S,W0,W1,Moment_Order);