function [y,z] = Get_HIST_General_Spatial(PARS,Model,Parameters,P0)
%% These are the five parameters for the Hog1 nuclearization curve.
Parameters = Define_Pars(Parameters,PARS,Model);

if nargin<4
    P0 = spalloc(PARS.N_states,1,1);
    P0(1)=1;
end
Jac = @(t,x)JAC_wHog(t,x,PARS,Parameters,Model);
Ode = @(t,x)ODE_wHog(t,x,PARS,Parameters,Model);

options = odeset('Jacobian',Jac,'AbsTol',1e-8,'JPattern',PARS.JPattern,'MaxStep',1);
% try
    [~,YOUT] = ode15s(Ode,PARS.Time_Array,P0,options);  %solve CME at given time points
% catch
%     save INTEG_ERR
%     YOUT = zeros(length(PARS.Time_Array),length(P0));
%     YOUT(:,1) = 1;
% end

if size(YOUT,1)<length(PARS.Time_Array)  %% If there was an integration error.
    save INTEG_ERR
    YOUT = zeros(length(PARS.Time_Array),length(P0));
    YOUT(:,1) = 1;
end

Dists_Mat = ones(PARS.N_Space(1),PARS.N_Space(2)+1,PARS.N_Space(3)+1);   % Map between index and state
y = zeros(length(PARS.Time_Array),PARS.N_Space(2)+1,PARS.N_Space(3)+1);
z = zeros(length(PARS.Time_Array),PARS.N_Space(1));
for i_time=1:length(PARS.Time_Array)
    Dists_Mat(:) = YOUT(i_time,:);
    y(i_time,:,:)=squeeze(sum(Dists_Mat,1));
    z(i_time,:) = squeeze(sum(sum(Dists_Mat,2),3));
end
y = max(y,1e-14);

function [Jac] = JAC_wHog(t_real,~,PARS,Parameters,Model)
t =max(0,t_real-Parameters.t_offset);
if t>0
    if isempty(Model.Signal.Function)
        hog_sig = (1-exp(-Parameters.r1*(t)))*exp(-Parameters.r2*(t));
        hog_sig = (hog_sig/(1+hog_sig*Parameters.M))^Parameters.eta/1.0727183838607e-10;
    else
        hog_sig = Model.Signal.Function(t_real);
    end
    
    if hog_sig>1e-8
        kprod = zeros(PARS.Num_Genes,PARS.Num_States);
        k=zeros(PARS.Num_States,PARS.Num_States);
        for i=1:PARS.Num_States
            for g=1:PARS.Num_Genes
                kprod(g,i) = max(-Model.Pars_Prods.Const(g,i),Model.Pars_Prods.Time(g,i)*hog_sig);
            end
            for j=1:PARS.Num_States
                if i~=j
                    k(i,j) =  max(-Model.Pars_State_Change.Const(i,j),Model.Pars_State_Change.Time(i,j)*hog_sig);
                end
            end
        end
        
        kdeg_nuc = max(-Model.Pars_Degrade_Nuc.Const(1),Model.Pars_Degrade_Nuc.Time(1)*hog_sig);
        
        Jac = Parameters.A_Constant+kdeg_nuc*PARS.Adeg_Nuc;
        for i=1:PARS.Num_States            
            for g=1:PARS.Num_Genes
                Jac = Jac + kprod(g,i)*PARS.Aprod(i,g).A;
            end
            for j=1:PARS.Num_States
                if i~=j
                    Jac = Jac + k(i,j)*PARS.A(i,j).A;
                end
            end
        end
                
        if strcmp(PARS.Pars_Type,'N_State_General')&&(strcmp(PARS.Spatial,'Spatial')||strcmp(PARS.Spatial,'Fixed'))
            kdeg_cyt = max(-Model.Pars_Degrade_Cyt.Const(1),Model.Pars_Degrade_Cyt.Time(1)*hog_sig);
            transp = max(-Model.Pars_Transport.Const(1),0);
            Jac = Jac+transp*PARS.ATransport + kdeg_cyt*PARS.Adeg_Cyt;
        elseif strcmp(PARS.Pars_Type,'Connected')&&strcmp(PARS.Spatial,'NonSpatial')&&PARS.Num_Genes==2
            kdeg_Sp2 = max(-Model.Pars_Degrade_Nuc.Const(2),Model.Pars_Degrade_Nuc.Time(2)*hog_sig);
            Jac = Jac+kdeg_Sp2*PARS.Adeg_Spec2;
        end
    else
        Jac = Parameters.A_Constant_NoSig;
    end
else
    Jac = Parameters.A_Constant_NoSig;
end
%
% if t>10
%     t
%     hog_sig
% [jacsum,ind] = max(abs(sum(Jac)))
% [Aconsum,ind] = max(abs(sum(Parameters.A_Constant)))
% Model.Pars_Prods.Const
% Model.Pars_Prods.Time
% pause
% end

function [fdot] = ODE_wHog(t_real,P,PARS,Parameters,Model)
fdot = JAC_wHog(t_real,[],PARS,Parameters,Model)*P;

function Parameters = Define_Pars(Parameters,PARS,Model)

Parameters.A_Constant = Model.Pars_Degrade_Nuc.Const(1)*PARS.Adeg_Nuc;
Parameters.A_Constant_NoSig = max(0,Model.Pars_Degrade_Nuc.Const(1))*PARS.Adeg_Nuc;
for i=1:PARS.Num_States
    for g=1:PARS.Num_Genes
        Parameters.A_Constant=Parameters.A_Constant +...
            Model.Pars_Prods.Const(g,i)*PARS.Aprod(i,g).A;
        Parameters.A_Constant_NoSig=Parameters.A_Constant_NoSig +...
            max(0,Model.Pars_Prods.Const(g,i))*PARS.Aprod(i,g).A;
    end
    for j=1:PARS.Num_States
        if i~=j;
            Parameters.A_Constant=Parameters.A_Constant +...
                Model.Pars_State_Change.Const(i,j)*PARS.A(i,j).A;
            Parameters.A_Constant_NoSig=Parameters.A_Constant_NoSig +...
                max(0,Model.Pars_State_Change.Const(i,j))*PARS.A(i,j).A;
        end
    end
end
if strcmp(PARS.Pars_Type,'N_State_General')&&(strcmp(PARS.Spatial,'Spatial')||strcmp(PARS.Spatial,'Fixed'))
    Parameters.A_Constant = Parameters.A_Constant+...
        Model.Pars_Transport.Const(1)*PARS.ATransport+...
        Model.Pars_Degrade_Cyt.Const(1)*PARS.Adeg_Cyt;
    Parameters.A_Constant_NoSig = Parameters.A_Constant_NoSig+...
        max(0,Model.Pars_Transport.Const(1))*PARS.ATransport+...
        max(0,Model.Pars_Degrade_Cyt.Const(1))*PARS.Adeg_Cyt;
elseif strcmp(PARS.Pars_Type,'Connected')&&strcmp(PARS.Spatial,'NonSpatial')&&PARS.Num_Genes==2
    %     error('Need to add second species degradation');
    Parameters.A_Constant = Parameters.A_Constant+...
        Model.Pars_Degrade_Nuc.Const(2)*PARS.Adeg_Spec2;
    Parameters.A_Constant_NoSig = Parameters.A_Constant_NoSig+...
        max(0,Model.Pars_Degrade_Nuc.Const(2))*PARS.Adeg_Spec2;
end


