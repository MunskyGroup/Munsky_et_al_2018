function make_mom_conv_fun(Ns,n,FUNNAME)
% Creates a code that converts from uncentered to centered moments.
%% Delete previous version if it exists.
if exist(FUNNAME,'file');
    disp('Moment conversion function already exists')
    
    B = zeros(1,Ns);
    nn = 0;
    for k = 1:(n+1)^Ns-1
        B(1) = B(1)+1;
        j=1;
        while B(j)==n+1
            B(j)=0;
            B(j+1) = B(j+1)+1;
            j=j+1;
        end
        if sum(B)<=n;
            nn=nn+1;
        end
    end
    if nargin(FUNNAME)==nn
%         disp('Previous file appears to be the right size')
    else
%         disp('Previous function appears to be incorrect')
        delete(FUNNAME)
        make_mom_conv_fun(Ns,n,FUNNAME)
    end
    
else
%     disp('Deriving moment conversion function')
    
    
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
    Bred =Bred';
    
    %% Create vector of all uncentered moments
    for k=2:K
        a = num2str(Bred(:,k)');
        a=['m',a(a~=' ')];
        uncent(k,1) = sym(a);
    end
    uncent(1) = 1;
    
    %% Find the centered moments in terms of the uncentered moments.
    for k = 2:K
        b = Bred(:,k);
        cent(k,1) = uncent(k);
        for ik=k-1:-1:1
            bdiff = b-Bred(:,ik);
            if min(bdiff)>=0
                pref=uncent(ik);
                for is=1:Ns
                    pref = pref*(-1)^bdiff(is)*...
                        nchoosek(b(is),bdiff(is))*...
                        uncent(is+1)^bdiff(is);
                end
                cent(k) = cent(k) + pref;
            end
        end
    end
    %% Save function to do this next time.
    cent = cent(2:end);
    matlabFunction(cent,'file',FUNNAME,'vars',uncent(2:end));
%     disp(['Successfully create function ',FUNNAME]);
%     disp(['for Ns=',num2str(Ns),' and n=',num2str(n)])
end