W0 = [2.3;0;0;0];
W1 = [0 0;1 0;1 0;0 1];
S = [1 -1 0  0;...
    0  0 1 -1];
n=2;

MOMFUN='mmfun.m';
if exist(MOMFUN,'file')
    delete(MOMFUN)
end

tic
if ~exist(MOMFUN,'file')
    get_mom_eqns_linear(W0,W1,S,n,MOMFUN);
end
toc

MOMFUN = str2func(MOMFUN(1:end-2));

X = [W0',W1(:)'];
for i=1:length(X); v{i} = X(i); end
[a,b] = MOMFUN(v{:})
