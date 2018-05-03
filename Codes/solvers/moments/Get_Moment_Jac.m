function Mom_Jacobian = Get_Moment_Jac(S,W0,W1,Moment_Order)
n_spec = size(W1,2);

%% First for the means S*W1*MU;  %Derivative of Means.
A_Mn  = S*W1;

if Moment_Order==1
    Mom_Jacobian = A_Mn;
elseif Moment_Order==2;
    A_Mn(n_spec+sum(1:n_spec),n_spec+sum(1:n_spec))=0;
    A_Mn = sparse(A_Mn);
    
    %% Now for the second moment
    % dSIG = S*W1*SIG + SIG*W1'*S' + S*diag((W1*MU + W0)')*S' + S*W0*MU'+MU*W0'*S';  % derivative of Co-Variances.
    %% S*W1*SIG and SIG*W1'*S'
    vec = zeros(n_spec,n_spec);
    vec(1,:) = ones(1,n_spec);
    for i=2:n_spec
        vec(i,1) = i;
        vec(i,2:end) = [(n_spec-1:-1:n_spec-i+1),ones(1,n_spec-i)];
    end
    
    vec_map = zeros(1,sum(1:n_spec));
    vec_map(1:n_spec)=(1:n_spec);
    for i=2:n_spec
        vec_map(1,sum(n_spec:-1:n_spec-(i-2))+1:sum(n_spec:-1:n_spec-(i-1))) = (i-1)*n_spec+(i:n_spec);
    end
    
    SW1 = S*W1;
    SW1_T = SW1';
    A_SW1 = zeros(n_spec+sum(1:n_spec));
    A_SW1_T = zeros(n_spec+sum(1:n_spec));
    itmp=0;
    for i=1:n_spec
        for j=i:n_spec
            itmp=itmp+1;
            A_SW1(  n_spec+itmp,n_spec+cumsum(vec(j,:))) = SW1(i,:);
            A_SW1_T(n_spec+itmp,n_spec+cumsum(vec(i,:))) = SW1_T(:,j);
        end
    end
    A_SW1=sparse(A_SW1);
    A_SW1_T=sparse(A_SW1_T);
    
    %% S*diag((W1*MU + W0)')*S'
    A_MU = zeros(n_spec+sum(1:n_spec));
    for i=1:n_spec
        TMP = S*diag(W1(:,i))*S';
        A_MU(n_spec+1:end,i) = TMP(vec_map);
    end
    A_MU=sparse(A_MU);
    
    %% S*W0*MU' and MU*W0'*S'
    A_SW0 = zeros(n_spec+sum(1:n_spec));
    TMP = S*W0*ones(1,n_spec);
    TMP_T = ones(n_spec,1)*W0'*S';
    itmp=n_spec;
    for j=1:n_spec
        for k=j:n_spec
            itmp=itmp+1;
            A_SW0(itmp,k)=TMP(j,k)+TMP_T(j,k);
        end
    end
    A_SW0 = sparse(A_SW0);
    
    %%
    Mom_Jacobian = (A_Mn + A_SW1 + A_SW1_T +  A_MU + A_SW0);
  
else
    n=Moment_Order;
    MOMFUN = ['mom_fun_order',num2str(n),'_spec',num2str(n_spec),'.m'];
    if ~exist(MOMFUN,'file')
        disp('Deriving moment equations -- this could take a while.')
        get_mom_eqns_linear(W0,W1,S,n,MOMFUN);
    end
    MOMFUN = str2func(MOMFUN(1:end-2));
    X = [W0',W1(:)'];
    for i=1:length(X); v{i} = X(i); end
    Mom_Jacobian = MOMFUN(v{:}); 
    Mom_Jacobian=Mom_Jacobian(2:end,2:end);
end