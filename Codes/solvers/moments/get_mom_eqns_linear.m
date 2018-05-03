function [A,B] = get_mom_eqns_linear(~,~,S,n,OLDFUN)
% n moment order.
% OLDFUN is the name of a function to be written in the form 
%     dMom/dt = fun(Mom);
%%
Ns = size(S,1);  % Number of species.
M = size(S,2);   % Number of reactions

w0 = sym('w0_',M,'positive');
w1 = sym('w1_',[M,N],'real');

%% Find all possible combinations of up to order n for each species.
BPos = zeros((n+1)^Ns,Ns);
B = zeros(1,Ns);
BPos(1,:) = B;
for k = 1:(n+1)^Ns-1
    B(1) = B(1)+1;
    j=1;
    while B(j)==n+1
        B(j)=0;
        B(j+1) = B(j+1)+1;
        j=j+1;
    end
    BPos(k+1,:) = B;
end

%% Reorder all possible combinations by order up to total order n.
K = sum(sum(BPos,2)<=n);
Bred = zeros(K,Ns);
Kmom=0;
nM = ones(n,1);
for i=0:n
    J = find(sum(BPos,2)==i);
    nM(i+1) = length(J);
    Bred(Kmom+1:Kmom+nM(i+1),:)=BPos(J,:);
    Kmom=Kmom+nM(i+1);
end

%% Find symbolic expressions for RHS of moment equations.
X = sym('x',[Ns,1],'integer');
for k=1:Kmom
    f = 0;
    for j=1:M
        w = w0(j) + w1(j,:)*X;
        s = S(:,j);
        f = expand(f+w*(prod((X+s).^(Bred(k,:)'))-prod(X.^(Bred(k,:)'))));
    end
    RHS(k,1) = f;
end

%% Convert to linear system dx/dt = A*x + B for UNCENTERED moments.
clear A B
for k=1:Kmom
    waitbar(k/Kmom)
    B(k,1) = subs(RHS(k,1),X,zeros(size(X)));
    for j=1:Kmom
        tmp = RHS(k,1);
        for i=1:Ns
            tmp = 1/factorial(Bred(j,i))*diff(tmp,Bred(j,i),X(i));
        end
        tmp = subs(tmp,X,zeros(size(X)));
        A(k,j) =tmp;
    end
end
matlabFunction(A,B,'file',OLDFUN,'vars',[w0,w1],'Sparse',true)