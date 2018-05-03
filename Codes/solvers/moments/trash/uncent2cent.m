function [cent] = uncent2cent(uncent,Ns,n)
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
Bred = Bred';

%% 
X = sym('x',[Ns,1],'integer');
for i=1:n
    mu(i).mu = sym(['mu',num2str(i)],[nM(i+1),1],'integer');
end
for k = 1:K
    b = Bred(:,k);
    nk = sum(b);
    if nk<2
        cent(k,1) = uncent(k);    
    else
        cent(k,1)=0;
        q=X-mu(1).mu;
        v = expand(prod(q.^b));
        for ik=k:-1:1     
            if v==0
                break
            end
            tmp=v;
            for i=1:Ns
                tmp = 1/factorial(Bred(i,ik))*diff(tmp,Bred(i,ik),X(i));
            end
            tmp = subs(tmp,X,zeros(size(X)));
            if tmp~=0
                cent(k) =  cent(k)+subs(tmp,mu(1).mu,uncent(2:Ns+1))*uncent(ik);
                v = v - tmp*prod(X.^Bred(:,ik));
            end
        end
    end
end
    
% X = sym('x',[Ns,1],'integer');
% mu1 = sym('x',[Ns,1],'integer');
% for k=1:Kmom
%     n_max = max(b);
%     f = expand(prod((X-).^(Bred(k,:)'));
%     
%     RHS(k,1) = f;
% end
% 


