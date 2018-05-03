function [Ts_all_MOD] = Model_TS_Results_SSA(SSA,MOD,speed,PD,N,f2,f3,f4,mod_col,NSamps,N_cells,inten_dist)
if nargin<10
    DownSamp = SSA.TS;
    N_cells = size(SSA.TS,1)*ones(1,16);
elseif nargin==10
    DownSamp = SSA.TS(ceil(size(SSA.TS,1)*rand(NSamps,1)),:,:);
    N_cells = NSamps*ones(1,16);
else
    DownSamp = SSA.TS(ceil(size(SSA.TS,1)*rand(NSamps,1)),:,:);
end
    
NRN = zeros(size(DownSamp,1),size(DownSamp,2));

% Find the mRNA who start at specific times in the past.
Nk = size(DownSamp,3);      % Number of SSA time points for tS.
dL = (MOD.gene_length/Nk)*(speed/MOD.k_elong);   % resolution of TS data.
psn = [Nk-1:-1:0]*dL;
Starts_TS = squeeze(DownSamp(:,:,2:Nk)-DownSamp(:,:,1:Nk-1));
CC = zeros(1,Nk);
for ik = 1:Nk
    if psn(ik)>MOD.gene_length
        CC(ik)=0;
    else
        CC(ik)=sum(PD<psn(ik));
    end
end
for ik = 1:Nk-1
    if CC(ik)>0
        NRN = NRN + squeeze(Starts_TS(:,:,ik))*CC(ik)/max(CC);
    end
end

%% Step 3b -- Plot the fraction of cells with TS versus time (MODEL).
frac_TS_v_t_MOD = sum(NRN>N)/size(DownSamp,1);
figure(f2);%subplot(3,1,1);
plot(MOD.tt,frac_TS_v_t_MOD,'r','linewidth',3); hold on;

%% Step 5b -- Plot the TS distributions versus time (MODEL).
figure(f3);
for it=1:length(MOD.tt)
    subplot(4,4,it);
    B = [N:40];
        TS_v = NRN(1:N_cells(it),it);
        [A,B] = hist(TS_v,B);
%         stairs(B,A/sum(A),'linewidth',3,'color',mod_col); hold on
    
%     A = inten_dist(1:length(B),it);
%     stairs(B,A/sum(A),'linewidth',3,'color',mod_col); hold on
    plot(B,A/sum(A),'r','linewidth',3); hold on
end

%% Step 6b -- Plot the TS distributions averaged over time (MODEL).
figure(f4);
B = [N:50];
A = zeros(length(B),1);
sm=0;
Ts_all_MOD=[];
for it=1:length(MOD.tt)
    TS_v = NRN(1:N_cells(it),it);
    Ts_all_MOD = [Ts_all_MOD,TS_v'];
     sm = sm+N_cells(it);
end
[A,B] = hist(Ts_all_MOD,B);
plot(B,A/sm,'r','linewidth',3); hold on

Ts_all_MOD = [N:49]*A(2:end)'/(sum(A(2:end)));

