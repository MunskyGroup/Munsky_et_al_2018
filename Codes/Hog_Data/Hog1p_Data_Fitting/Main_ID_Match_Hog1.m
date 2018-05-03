
function Main_ID_Match_Hog1
clear all
global n_iters min_err

load best_par_hog1_signal
x0=[Par(1:6)];
min_err=inf;
min_err = get_err(x0)
return
% x0=[2.3613e+00 1.2928e-05 1.6773e-03  1.5764e+04]

% x0=[2.3675 1.3045e-05 0.001677  15653  0.3 0.01824 53.249];

%% initialize Random number Generator
randn('state',sum(100*clock)); rand('twister',sum(100*clock));  %seed the randome number generator
n_iters=0;  %iteration counter
min_err = 80;  %initialize objective function

%%  Gradient search
options = optimset('MaxFunEvals', 40000, 'MaxIter', 10000,'FunValCheck', 'on');
for jj=1:10
    
    for i=1:2
        x0 = fminsearch(@get_err,x0,options);  %call gradient based search
        min_err
    end
    ann_gen = @(x,L) (x+x.*randn(size(x))/1000*min(L,1)^2);
    cooling = @(T) (0.5*T);
    options = anneal();
    options.Generator = ann_gen;
    options.CoolSched = cooling;
    options.Verbosity = 0;
    options.InitTemp = 1e-3;
    options.StopTemp = 1e-10;
    [x0,fval] = annealT(@get_err,x0,options);
end

function [err] = get_err(Par)
global min_err
Parameters.eta = Par(1);
Parameters.r1 = Par(2);
alpha = Par(3);
Parameters.M = Par(4);
Parameters.del = Par(5);
Parameters.t0 = Par(6);

load Hog1-02M-04M
Ht02=M02/max(M04);
Ht04=M04/max(M04);

t02 = max(0,t02*60-Parameters.t0);
t04 = max(0,t04*60-Parameters.t0);

Parameters.r2 = alpha/(.2-Parameters.del);
hog_sig = (1-exp(-Parameters.r1*t02)).*exp(-Parameters.r2*t02);
k_02 = (hog_sig./(1+hog_sig*Parameters.M)).^Parameters.eta/1.0727183838607e-10;

Parameters.r2 = alpha/(.4-Parameters.del);
hog_sig = (1-exp(-Parameters.r1*t04)).*exp(-Parameters.r2*t04);
k_04 = (hog_sig./(1+hog_sig*Parameters.M)).^Parameters.eta/1.0727183838607e-10;

err = norm(Ht02-k_02,2)^2+norm(Ht04-k_04,2)^2;

if err<min_err
%     Par
    save best_par_hog1_signal Par
    min_err=err;
    hold off
    plot(t02/60,k_02,'r',t04/60,k_04,'b','linewidth',3);
    hold on
    plot(t02/60,Ht02,'ro',t04/60,Ht04,'bo','linewidth',3);
    drawnow
    set(gca,'fontsize',16,'xlim',[-2 30],'ylim',[0 1.05])
    legend('0.2M, WT','0.4M, WT')
    
    save WT-02-04

end

