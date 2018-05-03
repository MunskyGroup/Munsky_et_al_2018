function Vals =cme_General_Spatial_NState(Vals)
%% This function will define the master equation of the stochastic process.

N = Vals.N_Space;  % States, Max Nuc, Max Cyt
N_total = Vals.N_states;  % Maximum number of states.

Indexing = ones(N(1),N(2)+1,N(3)+1);   % Map between index and state
Indexing(:) = (1:N_total);

for i=1:N(1)
    %% State Changes, Si -> Sj
    Ix = zeros(1*(N(2)+1)*(N(3)+1),1);
    Ix(:) = Indexing(i,:,:);  % Index of reaction start states. This reaction can occur when in state Si.
    for j=1:N(1)
        if i~=j
            Vals.A(i,j).A = sparse(Ix,      Ix,-1,N_total,N_total)+...     %
                sparse(Ix+(j-i),Ix,1 ,N_total,N_total);          % reaction changes the state by (j-i)
        end
    end
    %% Production of mRNA
    Ixp = zeros(N(2)*(N(3)+1),1);
    Ixp(:) = Indexing(i,1:N(2),:);  % This prod reaction can occur when in state Si, and can result in up to N(2) mRNA's of species 1.
    Vals.Aprod(i,1).A = sparse(Ixp,     Ixp,-1,N_total,N_total)+...
                        sparse(Ixp+N(1),Ixp, 1,N_total,N_total); % reaction changes the state by +(N(1))   
end

%% Degradation of mRNA (Nuclear or Species 1).
Ix = zeros(N(1)*N(2)*(N(3)+1),1);
Ix(:)= Indexing(:,2:N(2)+1,:);  % Index of reaction start states. This reaction can occur for all states except when there are no nuclear mRNA
[~,i2,~] = ind2sub([N(1),N(2)+1,N(3)+1],Ix);   % i2-1 is the number of nuclear mRNA.
i2 = i2-1;  % change to number of nuclear mRNA
Vals.Adeg_Nuc = sparse(Ix,     Ix,-i2,N_total,N_total)+...       % reaction rate scales with the number of nuclear mRNA.
                sparse(Ix-N(1),Ix,i2,N_total,N_total); % reaction changes the state by -N(1)

Vals.JPattern = sparse(abs(Vals.Adeg_Nuc));
                      
if strcmp(Vals.Pars_Type,'N_State_General')&&strcmp(Vals.Spatial,'NonSpatial')
    Vals.Adeg_Cyt=0*Vals.JPattern;
    Vals.ATransport=Vals.JPattern;
elseif strcmp(Vals.Pars_Type,'N_State_General')&&(strcmp(Vals.Spatial,'Spatial')||strcmp(Vals.Spatial,'Fixed'))
    %% Transport of Nuclear mRNA to Cytoplasm.
    Ix = zeros(N(1)*N(2)*N(3),1);
    Ix(:)= Indexing(:,2:N(2)+1,1:N(3));  % Index of reaction start states. This reaction can occur for all states except when there are no nuclear mRNA - or - when there are maximal cyt mRNA.
    [~,i2,~] = ind2sub([N(1),N(2)+1,N(3)+1],Ix);   % i2-1 is the number of nuclear mRNA.
    i2 = i2-1;  % change to number of nuclear mRNA
    Vals.ATransport = sparse(Ix,Ix,-i2,N_total,N_total)+...       % reaction rate scales with the number of nuclear mRNA.
        sparse(Ix+N(1)*N(2),Ix,i2,N_total,N_total); % reaction changes the state by -N(1)+N(1)*(N(2)+1) = N(1)*N(2)
    Vals.JPattern = Vals.JPattern+abs(Vals.ATransport);
    %% Degradation of Cytoplasmic mRNA.
    Ix = zeros(N(1)*(N(2)+1)*N(3),1);
    Ix(:)= Indexing(:,:,2:N(3)+1);  % Index of reaction start states. This reaction can occur for all states except when there are no cytoplasmic mRNA
    [~,~,i3] = ind2sub([N(1),N(2)+1,N(3)+1],Ix);   % i3-1 is the number of cyt. mRNA.
    i3 = i3-1;  % change to number of cyt. mRNA
    Vals.Adeg_Cyt = sparse(Ix,              Ix,-i3,N_total,N_total)+...       % reaction rate scales with the number of cyt. mRNA.
                    sparse(Ix-N(1)*(N(2)+1),Ix, i3,N_total,N_total); % reaction changes the state by -N(1)*(N(2)+1)
    Vals.JPattern = Vals.JPattern+abs(Vals.Adeg_Cyt);

elseif strcmp(Vals.Pars_Type,'Connected')&&strcmp(Vals.Spatial,'NonSpatial')&&Vals.Num_Genes==2
    % We now need to add in the events for the second gene's transcription
    % and degradation events.
    for i=1:N(1)
        %% Production of mRNA, species 2
        Ixp = zeros((N(2)+1)*(N(3)),1);
        Ixp(:) = Indexing(i,:,1:N(3));  % This prod reaction can occur when in state Si, and can result in up to N(3) mRNA's of species 2.
        Vals.Aprod(i,2).A = sparse(Ixp,     Ixp,-1,N_total,N_total)+...
                            sparse(Ixp+N(1)*(N(2)+1),Ixp, 1,N_total,N_total); % reaction changes the state by +(N(1)*(N(2)+1))
    end
    Ix = zeros(N(1)*(N(2)+1)*N(3),1);
    Ix(:)= Indexing(:,:,2:N(3)+1);  % Index of reaction start states. This reaction can occur for all states except when there are no cytoplasmic mRNA
    [~,~,i3] = ind2sub([N(1),N(2)+1,N(3)+1],Ix);   % i3-1 is the number of cyt. mRNA.
    i3 = i3-1;  % change to number of cyt. mRNA
    Vals.Adeg_Spec2 = sparse(Ix,              Ix,-i3,N_total,N_total)+...       % reaction rate scales with the number of cyt. mRNA.
                      sparse(Ix-N(1)*(N(2)+1),Ix, i3,N_total,N_total); % reaction changes the state by -N(1)*(N(2)+1)
else
    Vals
    error('FSP not set up for this situation')
end

for i=1:N(1)
    for g=1:Vals.Num_Genes
        Vals.JPattern = Vals.JPattern+sparse(abs(Vals.Aprod(i,g).A));
    end
    for j=1:N(1)
        if i~=j;
            Vals.JPattern = Vals.JPattern+sparse(abs(Vals.A(i,j).A));
        end
    end
end
Vals.JPattern = sparse(Vals.JPattern~=0);
