function PARS = PARS_Definition_Connected_Models(PARS)
PARS.Prod = zeros(PARS.Num_Genes,PARS.Num_States);
PARS.Prod_Time= zeros(PARS.Num_Genes,PARS.Num_States);
PARS.Degrade = zeros(PARS.Num_Genes,2);
PARS.Degrade_Time = zeros(PARS.Num_Genes,2);
PARS.Transport = zeros(PARS.Num_Genes,1);
PARS.Transport_Time = zeros(PARS.Num_Genes,1);
PARS.States = zeros(PARS.Num_States,PARS.Num_States);
PARS.States_Time = zeros(PARS.Num_States,PARS.Num_States);

switch PARS.Model_Connection
    case 0
    case '4State'  % 4-state connected model
        PARS.Num_States = 4;
        Trans_States = [1 2 3 4;...
                        1 2 3 4];
    case '6StateLow' % 6-state-Low connected model
        PARS.Num_States = 6;
        Trans_States = [1 2 1 2 3 4;...
                        1 1 2 2 3 4];
    case '6StateMid' % 6-state-Mid connected model
        PARS.Num_States = 6;
        Trans_States = [1 2 3 2 3 4;...
                        1 2 2 3 3 4];
    case '6StateHigh' % 6-state-High connected model
        PARS.Num_States = 6;
        Trans_States = [1 2 3 4 3 4;...
                        1 2 3 3 4 4];
    case '8State'
        PARS.Num_States = 8;
        Trans_States = [1 2 1 2 3 4 3 4;...
                        1 1 2 2 3 3 4 4];
    case '10StateLow'
        PARS.Num_States = 10;
        Trans_States = [1 2 3 1 2 3 1 2 3 4;...
                        1 1 1 2 2 2 3 3 3 4];
    case '10StateHigh'
        PARS.Num_States = 10;
        Trans_States = [1 2 3 4 2 3 4 2 3 4;...
                        1 2 2 2 3 3 3 4 4 4];
end

ipar = 0;
PARS.Spatial_Pars = [];
% State Transitions
for j = 1:PARS.Num_States
    for k = 1:PARS.Num_States
        if j~=k
            ipar = ipar+1;
            PARS.States(j,k) = PARS.Vector(ipar);
        end
    end
end
% State Transitions in Time
for j = 1:PARS.Num_States
    for k = 1:PARS.Num_States
        if j~=k
            ipar = ipar+1;
            PARS.States_Time(j,k) = PARS.Vector(ipar);
        end
    end
end
% Production
for i = 1:PARS.Num_Genes
    for j = 1:max(max(Trans_States))
        ipar = ipar+1;
        PARS.Prod(i,(Trans_States(i,:)==j)) = PARS.Vector(ipar);
    end
    % Production Time
    for j = 1:max(max(Trans_States))
        ipar = ipar+1;
        PARS.Prod_Time(i,(Trans_States(i,:)==j)) = PARS.Vector(ipar);
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