function [speed] = get_TS_Fit_errs_All(RNAdata,gene,inacl,reps,spat,type)
N = 2; % Minimum number of nascent mRNA required to be considered as a TS.

if strcmp(gene,'BOTH')
    [DAT1,~,TS_ALL1,~] = get_DAT(RNAdata,'STL1',N,inacl);
    [DAT2,~,TS_ALL2,~] = get_DAT(RNAdata,'CTT1',N,inacl);
    gt_err = @(x)(get_Fit_errs_FSP_TS(10.^x,inacl,reps,'STL1',TS_ALL1,N,DAT1,spat,type)+...
        get_Fit_errs_FSP_TS(10.^x,inacl,reps,'CTT1',TS_ALL2,N,DAT2,spat,type));
else
    [DAT,~,TS_ALL,~] = get_DAT(RNAdata,gene,N,inacl);
    gt_err = @(x)get_Fit_errs_FSP_TS(10.^x,inacl,reps,gene,TS_ALL,N,DAT,spat,type);
end
options = optimset('display','iter');
speed = 10.^(fminsearch(gt_err,2,options));

% speed_ar = logspace(log10(30),log10(200),100);
% obj = 0*speed_ar;
% for i=1:length(speed_ar)
%     obj(i) = gt_err(log10(speed_ar(i)));
%     if strcmp(gene,'STL1')
%         figure(101);clf;
%     else
%         figure(102);clf;
%     end
%     loglog(speed_ar,obj,'x');
%     drawnow
% end
% [~,j]=min(obj);
% speed = speed_ar(j);

function ERR = get_Fit_errs(speed,DAT,MOD,Starts,N,TS_ALL,on_vs_t_dat,inacl,gene)
load('Hog_Data/Hog_Probes.mat')

switch gene
    case 'STL1'
        PD = PD_STL1;
    case 'CTT1'
        PD = PD_CTT1;
end

Nk = size(Starts(inacl(1)).TS,3);      % Number of SSA time points for tS.
ERR = 0;

TS_ALL_Mod = [];
on=zeros(2,16);
off=zeros(2,16);
for i_nacl = inacl
    %% Step 2b -- Find the TS spots intensities (MODEL)
    switch i_nacl
        case 1
            n_reps = 2;
        case 2
            n_reps = 3;
    end
    NRN = zeros(size(Starts(i_nacl).TS,1),size(Starts(i_nacl).TS,2));    
    % Find the mRNA who start at specific times in the past.
    dL = (MOD.gene_length/Nk)*(speed/MOD.k_elong);   % resolution of TS data.
    psn = [Nk-1:-1:0]*dL;
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
            NRN = NRN + squeeze(Starts(i_nacl).TS(:,:,ik))*CC(ik)/max(CC);
        end
    end
    %% Step 5b -- Plot the TS distributions versus time (MODEL).
    for it=1:length(MOD.tt)
        TS_v = NRN(:,it);
        [A,B] = hist(TS_v,[N:100]);
        A=A/sum(A);
        model(i_nacl,it).A = max(A,1e-5);
%         mod_mn(it) = mean(TS_v(TS_v>=N));        
        TS_ALL_Mod = [TS_ALL_Mod;TS_v];
            on(i_nacl,it)  = on(i_nacl,it)+sum(TS_v>=N);
            off(i_nacl,it) = off(i_nacl,it)+sum(TS_v<N);
    end
    
    %%
%         for j=1:n_reps
%             for it=2:16
%                 for k = 1:length(B)
%                     if ~isempty(DAT(i_nacl,j,it).A)&&DAT(i_nacl,j,it).A(k)>0
%                         ERR = ERR-DAT(i_nacl,j,it).A(k)*log(model(i_nacl,it).A(k));
%                     end
%                 end
%             end
%         end
    
%         for j=1:n_reps
%             for it=2:16
%                 if ~isempty(DAT(i_nacl,j,it).A)
%                     %                 whos
%                     %                 size(DAT(i_nacl,j,it).A)
%                     %                 size([N:length(DAT(i_nacl,j,it).A)])
%                     dat_mn = DAT(i_nacl,j,it).A(B>=N)*B(B>=N)';
% %                     [dat_mn,mod_mn(it)]
%                     %                  pause
%                     if ~isnan(mod_mn(it))
%                         ERR =ERR+(mod_mn(it)-dat_mn)^2;
%                     end
%                     
% %                     if isnan(ERR)
% %                         save TMP
% %                         return
% %                     end
%                 end
%             end
%         end
        
end
dist_TS_ALL = hist(TS_ALL,[0:100]);
dist_TS_ALL_Mod = hist(TS_ALL_Mod,[0:100]);
dist_TS_ALL_Mod=dist_TS_ALL_Mod/sum(dist_TS_ALL_Mod);
dist_TS_ALL_Mod=max(dist_TS_ALL_Mod,1e-6);
ERR = -sum(dist_TS_ALL.*log(dist_TS_ALL_Mod));
% 
% dist_TS_ALL = hist(TS_ALL,[N+1:50]);
% dist_TS_ALL_Mod = hist(TS_ALL_Mod,[N+1:50]);
% dist_TS_ALL_Mod=dist_TS_ALL_Mod/sum(dist_TS_ALL_Mod);
% dist_TS_ALL_Mod=max(dist_TS_ALL_Mod,1e-6);
% ERR = -sum(dist_TS_ALL.*log(dist_TS_ALL_Mod));

% on_vs_t_mod = on./(on+off);
% ERR = sum(sum((on_vs_t_mod(inacl,:)-on_vs_t_dat(inacl,:)).^2));
 
% mn_on_mod = mean(TS_ALL_Mod(TS_ALL_Mod>N));
% mn_on_dat = mean(TS_ALL(TS_ALL>N));
% ERR = (mn_on_mod-mn_on_dat)^2;
