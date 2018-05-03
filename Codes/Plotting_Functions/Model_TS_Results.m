function [Ts_all_MOD] = Model_TS_Results(SSA,MOD,speed,PD,N,f2,f3,f4,mod_col,NSamps,N_cells,inten_dist,fnas)
%% Step 3b -- Plot the fraction of cells with TS versus time (MODEL).
figure(f2);
frac_TS_v_t_MOD_FSP = 1-inten_dist(1,:);
plot(MOD.tt,frac_TS_v_t_MOD_FSP,'linewidth',3,'color',[0 0 1]); hold on;
%% Step 5b -- Plot the TS distributions versus time (MODEL).
figure(f3);
for it=1:length(MOD.tt)
    subplot(4,4,it);
    B = [N:40];
    A = inten_dist(1:length(B),it);
    plot(B,A,'linewidth',3,'color',mod_col); hold on
    set(gca,'ylim',[1e-4 1],'yscale','log','xlim',[0 25],'yticklabel',[])
end

%% Step 5b -- Plot the TS distributions versus time (MODEL).
figure(fnas);
% tar = [1 5 7 8 9 11 13]; 
tar = [5 6 7 9 11]; 
for it=1:length(tar)
    subplot(1,length(tar),it);
    B = [N:40];
    A = inten_dist(1:length(B),tar(it));
    plot(B,A,'linewidth',3,'color',mod_col); hold on
    set(gca,'ylim',[1e-4 1],'yscale','log','xlim',[0 25],'yticklabel',[])
    if it==1
        set(gca,'ytick',10.^([-4 -2 0]))
    end
end

%% Step 6b -- Plot the TS distributions averaged over time (MODEL).
figure(f4);
B = [N:50];
A = zeros(length(B),1);
sm=0;
for it=1:length(MOD.tt)
    A = A + inten_dist(1:length(B),it)*N_cells(it);
    sm = sm+N_cells(it);
end
plot(B,A/sm,'linewidth',3,'color',mod_col); hold on

Ts_all_MOD = [N:49]*A(2:end)/(sum(A(2:end)));

