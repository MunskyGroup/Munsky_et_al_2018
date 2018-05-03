function Reactions = assign_reactions_nstate_class(Reactions,obj,PARS) 
% Reactions - A Structure that will be used to contain all of the different
% reaction types.  The substructures of this will be:
%   Pars_State_Change   - Transitions from one state to another.
%   Pars_Prods          - Productionn of nuclear mRNA.
%   Pars_Degrade_Nuc    - Degradation of nuclear mRNA.
%   Pars_Degrade_Cyt    - Degradation of cytoplasmic mRNA.
%   Pars_Transport      - Degradation of nuclear to cytoplasmic mRNA.
%
% The input PARS contains the parameters for the current set of reactions.
% These are listed in the substructure 'States', which is a
% 3-dimensional matrix (Mgenes,Nstates,Nstates).  The first dimension corresponds to the particular
% gene number.  The second and third dimensions are the state or
% destination and the state of origination for each reaction. These will be
% combined into a large matrix of size N^M by N^M.

Reactions.Pars_State_Change.Const=zeros(obj.Num_States^obj.MGenes,obj.Num_States^obj.MGenes);
Reactions.Pars_State_Change.Time=zeros(obj.Num_States^obj.MGenes,obj.Num_States^obj.MGenes);

RXNs = zeros(obj.Num_States^obj.MGenes*obj.Num_States*obj.MGenes,3);
RXNs_Time = zeros(obj.Num_States^obj.MGenes*obj.Num_States*obj.MGenes,3);
k_ind = 0;

for k = 1:obj.Num_States^obj.MGenes
    ST_cell = cell([1,obj.MGenes]);
    [ST_cell{:}] = ind2sub(obj.Num_States*ones(1,obj.MGenes),k);
    ST = [ST_cell{:}]; % Here 'ST' is the state of each gene. 
    % We wish to find the neighbors of this particular point and the rates
    % at which the cells should transition to these neighbors.
    for i=1:obj.MGenes
        Orig = ST(i);    % The current state of gene i.
        for j = 1:obj.Num_States % The next state for gene i.
            if j~=Orig; % If there is a difference.
                k_ind=k_ind+1; % Increment reaction counter.
                RXNs(k_ind,:) = [k,k+(j-Orig)*obj.Num_States^(i-1),PARS.States(i,Orig,j)];
                % Record this reaction.
                RXNs_Time(k_ind,:) = [k,k+(j-Orig)*obj.Num_States^(i-1),PARS.States_Time(i,Orig,j)];
            end
        end  
        Reactions.Pars_Prods.Const(i,k)=PARS.Prod(i,Orig); 
        Reactions.Pars_Prods.Time(i,k)=PARS.Prod_Time(i,Orig);
    end
end
RXNs = RXNs(1:k_ind,:);
RXNs_Time = RXNs_Time(1:k_ind,:);

% Now to convert all of the reactions into a form that can be used later in
% the program.
for i=1:size(RXNs,1)
    Reactions.Pars_State_Change.Const(RXNs(i,1),RXNs(i,2))=RXNs(i,3);
    Reactions.Pars_State_Change.Time(RXNs_Time(i,1),RXNs_Time(i,2))=RXNs_Time(i,3);
end

for i=1:obj.MGenes
    Reactions.Pars_Degrade_Nuc.Const(i)=PARS.Degrade(i,1);
    Reactions.Pars_Degrade_Nuc.Time(i)=PARS.Degrade_Time(i,1);
    
    Reactions.Pars_Degrade_Cyt.Const(i)=PARS.Degrade(i,2);
    Reactions.Pars_Degrade_Cyt.Time(i)=PARS.Degrade_Time(i,2);
    
    Reactions.Pars_Transport.Const(i)=PARS.Transport(i);
    Reactions.Pars_Transport.Time(i)=PARS.Transport_Time(i);
end