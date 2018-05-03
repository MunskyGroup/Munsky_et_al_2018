function PARS = Reformat_PARS_Definition_Vector_Class(PARS)
PARS.Prod = zeros(PARS.Num_Genes,PARS.Num_States);
PARS.Prod_Time= zeros(PARS.Num_Genes,PARS.Num_States);
PARS.Degrade = zeros(PARS.Num_Genes,2);
PARS.Degrade_Time = zeros(PARS.Num_Genes,2);
PARS.Transport = zeros(PARS.Num_Genes,1);
PARS.Transport_Time = zeros(PARS.Num_Genes,1);
PARS.States = zeros(PARS.Num_Genes,PARS.Num_States,PARS.Num_States);
PARS.States_Time = zeros(PARS.Num_Genes,PARS.Num_States,PARS.Num_States);

ipar = 0;
PARS.Spatial_Pars = [];
for i=1:PARS.Num_Genes
    % State Transitions
    for j = 1:PARS.Num_States
        for k = 1:PARS.Num_States
            if j~=k
                ipar = ipar+1;
%                 PARS.States(i,j,k) = max(0,PARS.Vector(ipar));
                PARS.States(i,j,k) = PARS.Vector(ipar);
            end
        end
    end
    % State Transitions in Time
    for j = 1:PARS.Num_States
        for k = 1:PARS.Num_States
            if j~=k
                ipar = ipar+1;
                PARS.States_Time(i,j,k) = PARS.Vector(ipar);
            end
        end
    end
    % Production
    for j = 1:PARS.Num_States
        ipar = ipar+1;
%         PARS.Prod(i,j) = max(0,PARS.Vector(ipar));
        PARS.Prod(i,j) = PARS.Vector(ipar);
    end
    % Production Time
    for j = 1:PARS.Num_States
        ipar = ipar+1;
        PARS.Prod_Time(i,j) = PARS.Vector(ipar);
    end
    
    % Degrade Nuc
    ipar = ipar+1;
%     PARS.Degrade(i,1) = max(0,PARS.Vector(ipar));
    PARS.Degrade(i,1) = PARS.Vector(ipar);
    % Degrade Nuc Time
    ipar = ipar+1;
    PARS.Degrade_Time(i,1) = PARS.Vector(ipar);
    
    % Degrade Cyt
    ipar = ipar+1;
%     PARS.Degrade(i,2) = max(0,PARS.Vector(ipar));
    PARS.Degrade(i,2) = PARS.Vector(ipar);
    PARS.Spatial_Pars = [PARS.Spatial_Pars,ipar];

    % Degrade Cyt Time
    ipar = ipar+1;
    PARS.Degrade_Time(i,2) = PARS.Vector(ipar);
    PARS.Spatial_Pars = [PARS.Spatial_Pars,ipar];
    
    % Transport 
    ipar = ipar+1;
%     PARS.Transport(i) = max(0,PARS.Vector(ipar));
    PARS.Transport(i) = PARS.Vector(ipar);
    PARS.Spatial_Pars = [PARS.Spatial_Pars,ipar];

    % Transport 
    ipar = ipar+1;
    PARS.Transport_Time(i) = PARS.Vector(ipar);
    PARS.Spatial_Pars = [PARS.Spatial_Pars,ipar];

end
% PARS.Num_Genes
% PARS.States
% PARS.States_Time
% PARS.Prod
% PARS.Prod_Time
% PARS.Degrade
% PARS.Degrade_Time
% PARS.Transport
% PARS.Transport_Time
% PARS.Vector
%  pause

% Offset Time.
ipar = ipar+1;
PARS.t_offset = PARS.Vector(ipar);