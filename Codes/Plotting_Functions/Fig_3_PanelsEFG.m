function Fig_3_PanelsEFG(Gene)
ylim = [1e-4 1e2];
load(['Fit_Results/',Gene,'_0p2and0p4/NonSpatial/Distributions/4State_k21only_noloop/Analyses.mat']);
whos

W=0.13;
N = length(Hog.Distributions(1).distributions)-1;
Nt=12;

figure(11);
set(gcf,'Position',[1000 1172 769 166]);
subplot(1,5,1);
set(gca,'position',[0.070 0.1313 0.15 0.7310])

P = Hog.Distributions(1).distributions(Nt,:);
mu = Hog.Distributions(1).distributions(Nt,:)*[0:N]';
mu2 = Hog.Distributions(1).distributions(Nt,:)*[0:N]'.^2;
plot([0:N],P,'k','linewidth',2); hold on;
mu_old = mu; mu_new = mu;

F = cumsum(P);F=F/F(end);
I = find(F<0.5,1,'last');

if isempty(I); I =1;end
% plot((I-1)*[1 1],[1e-20,100],'b--')
sig = sqrt(mu2-mu^2);
vv = linspace(-N/2,N,1000);
stp = vv(2)-vv(1);

plot(vv,normpdf(vv,mu,sig)*stp,'m','linewidth',2);hold on
set(gca,'fontsize',14,'yscale','log','ylim',ylim,'ytick',10.^(-4:2:2),'xlim',[-10 40],'xtick',[0 15 30] );
plot([mu mu],[1e-20,100],'m--')
% plot((I-1)*[1 1],[1e-20,100],'b--')
title('1 cell');
% ylabel('Probability Density')

figure(12);
set(gcf,'Position',[1000         932         769         166]);
subplot(1,5,1)
set(gca,'position',[0.070    0.1313    W   0.7310])
plot([0:N],P,'k','linewidth',2); hold on;

plot(vv,normpdf(vv,mu,sig)*stp,'m','linewidth',2);hold on;
plot(vv,normpdf(vv,mu,sig)*stp,'c--','linewidth',2);hold on;
set(gca,'fontsize',14,'yscale','lin','xlim',[-10 10]);
plot([mu mu],[1e-20,max(max(P),max(normpdf(vv,mu,sig)))],'m--')
title('1 cell');
% ylabel('Probability Density')

%%
for kl = 2:5
    switch kl
        case 2
            Nsum=100;
            xlim = [-1 2];
            xtick = [0 2 4];
            xlimln = [-1.5 4.5];
            xtickln = [-1,0,1,2];
            plsp ='-';
        case 3
            Nsum=300;
            xlim = [-0.2 1];
            xtick = [0 1 2];
            xlimln = [-0.5 2.25];
            xtickln = [0 0.5 1];
            plsp ='-';
        case 4
            Nsum=1000;
            xlim = [0 0.75];
            xtick = [-0,0.5,1.0];
            xlimln = [-0.25 1.25];
            xtickln = [0,0.25,0.5];
            plsp = '-';
        case 5
            Nsum=3000;
            xlim = [0 0.5];
            xtick = [0 0.25 0.5];
            xlimln = [-0.1 0.75];
            xtickln = [0,0.25,0.5];
            plsp = '-';
        case 6
            Nsum=10000;
            %                 xlim = [-2 6]*mu;
            xtick = [0,0.25,0.5,0.75];
            xlimln = [0 0.5];
            xtickln = [0,0.25,0.5];
            plsp = '-';
    end
    
    xlimplt = [floor(mu-100*sig/sqrt(Nsum)) ceil(mu+100*sig/sqrt(Nsum))];
    
    figure(11);
    subplot(1,5,kl);
    set(gca,'position',[0.070+0.18*(kl-1)    0.1313    W   0.7310])
    
    X = (0:N);
    mx= N;
    pad_xx = 0;
    pad_xx_sc = 0;
    
    if~isempty(P)
        
        PX = P;
        PXi = PX;
        for i=1:Nsum-1
            PXi = conv(PX,PXi);
            kmin = find(PXi/max(PXi)>1e-30,1,'first');
            kmax = find(PXi/max(PXi)>1e-30,1,'last');
            PXi = PXi(kmin:kmax);
            pad_xx = pad_xx+(X(2)-X(1))*(kmin-1);
            
        end
        X_sum = linspace(pad_xx,pad_xx+(X(2)-X(1))*(length(PXi)-1),length(PXi))/Nsum;
        PX_sum = PXi/(X_sum(2)-X_sum(1));
        
        
        plot(X_sum,PX_sum,'k','linewidth',2);hold on;
    end
    sig = sqrt((mu2-mu^2))/sqrt(Nsum);
    vv = linspace(xlimplt(1),xlimplt(2),1000);stp = vv(2)-vv(1);
    plot(vv,normpdf(vv,mu,sig),'m','linewidth',2);hold on;
    set(gca,'fontsize',14,'yscale','log','ylim',ylim,'ytick',10.^(-4:2:2),'xlim',xlimln,'xtick',xtick,'yticklabels',[]);
    plot([mu mu],[1e-20,100],'m--')
    F = cumsum(PX_sum);F=F/F(end);
    I = find(F<0.5,1,'last');
    plot(X_sum(I)*[1 1],[1e-20,max(max(PX_sum),max(normpdf(vv,mu,sig)))],'b--')
    title([num2str(Nsum),' cells']);
    
    figure(12);
    subplot(1,5,kl);
    set(gca,'position',[0.070+0.18*(kl-1)    0.1313    W   0.7310])
    if~isempty(P)
        plot(X_sum,PX_sum,'k','linewidth',2);hold on;
    end
    plot(vv,normpdf(vv,mu,sig),'m','linewidth',2);hold on;
    plot(vv,normpdf(vv,mu,sig),'c--','linewidth',2);hold on;
    set(gca,'fontsize',14,'yscale','lin','xlim',xlimplt,'xtick',xtickln);
    set(gca,'fontsize',14,'yscale','lin','xlim',xlim);
    if~isempty(P)
        plot([mu mu],[1e-20,max(max(PX_sum),max(normpdf(vv,mu,sig)))],'m--')
        F = cumsum(PX_sum);F=F/F(end);
        I = find(F<0.5,1,'last');
%         plot(X_sum(I)*[1 1],[1e-20,max(max(PX_sum),max(normpdf(vv,mu,sig)))],'b--')
    end
    title([num2str(Nsum),' cells']);
    
    drawnow
end
%%
% figure(13);clf
% set(gcf,'Position',[1000 1172 769 166]);
% Nt=12;
% P = Hog.Distributions(1).distributions(Nt,:);
% % P = normpdf([0:300],40,5);
% % NsumV = [1,100,300,1000,3000];
% NsumV = [1,1000,10000,100000];
% for kl=2:4
%     Nsum = NsumV(kl);
%     mn = P*([0:N])';
%     mn2 = P*([0:N].^2)';
%     sig2true = mn2-mn^2;
%     
%     Ns = 10000;
%     var_s = zeros(1,Ns);
%     for i=1:Ns
%         histo = sample_d(P,Nsum)/Nsum;
%         mn = histo*([0:N])';
%         mn2 = histo*([0:N].^2)';
%         var_s(i) = (mn2-mn^2);
%     end
%     
%     set(gcf,'position',[ 1000         692         769         166]);
%     subplot(1,5,kl);
%     set(gca,'position',[0.070+0.18*(kl-1)    0.1313    W   0.7310])
%     
%     [A,B] = hist(var_s,linspace(0,21,50));
% %     [A,B] = hist(var_s,50);
%     plot(B,A/sum(A),'k','linewidth',2);hold on;
%     plot(sig2true*[1 1],[0 0.1],'m--')
%     set(gca,'xlim',[0 20],'ylim',[0 0.1]);
%     set(gca,'fontsize',14)
%     mn = P*([0:N])';
%     mn2 = P*([0:N].^2)';
%     sig2true = mn2-mn^2;
%     
%     g = (sqrt(Nsum));
%     Xi2_range = (Nsum/g-1)*B/(sig2true/g);
%     Q = chi2pdf(Xi2_range/g,Nsum/g-1);
%     plot(B,Q/sum(Q),'c','linewidth',2)
%     
%     s1111 = P*(([0:N]-mn).^4)';
%     mu4 = s1111-(Nsum-3)/(Nsum-1)*sig2true^2;
%     Q = normpdf(B,sig2true,sqrt(1/Nsum*mu4));
%     plot(B,Q/sum(Q),'m','linewidth',2)
%     title([num2str(Nsum),' cells']);
% end
%%
for Nt = 1:16
    P = Hog.Distributions(1).distributions(Nt,:);
    mn = P*([0:N])';
    sig2true = P*(([0:N]-mn).^2)';
    mn4 = P*(([0:N]-mn).^4)';
    var_var_1 = @(NN)(sqrt(1/NN*(mn4-(NN-3)/(NN-1)*sig2true^2))-sig2true/100)^2;
    var_var_10 = @(NN)(sqrt(1/NN*(mn4-(NN-3)/(NN-1)*sig2true^2))-sig2true/10)^2;
    var_dev_1 = @(NN)(sqrt(sqrt(1/NN*(mn4-(NN-3)/(NN-1)*sig2true^2)))-sqrt(sig2true)/100)^2;
    var_dev_10 = @(NN)(sqrt(sqrt(1/NN*(mn4-(NN-3)/(NN-1)*sig2true^2)))-sqrt(sig2true)/10)^2;
    var_mn_1 = @(NN)(sqrt(sig2true/(NN-1))-mn/100)^2;
    var_mn_10 = @(NN)(sqrt(sig2true/(NN-1))-mn/10)^2;
    Nv_1(Nt) = fminsearch(var_var_1,100);
    Nv_10(Nt) = fminsearch(var_var_10,100);
    Nd_1(Nt) = fminsearch(var_dev_1,100);
    Nd_10(Nt) = fminsearch(var_dev_10,100);
    Nm_1(Nt) = fminsearch(var_mn_1,100);
    Nm_10(Nt) = fminsearch(var_mn_10,100);
%     var_var = @(NN)(sqrt(1/NN*(mn4-(NN-3)/(NN-1)*sig2true^2))-sig2true/100)^2;
%     Np2(Nt) = fminsearch(var_var,100);
end
figure(14);clf;
subplot(2,1,2)
semilogy([1 2 4 6 8 10 15 20 25 30 35 40 45 50 55],Nv_1(2:end),'k','linewidth',3); hold on
semilogy([1 2 4 6 8 10 15 20 25 30 35 40 45 50 55],Nv_10(2:end),'r','linewidth',3); hold on
set(gca,'ylim',[1e2,1e8],'fontsize',14,'ytick',10.^(2:2:8),'xlim',[0 55],'xtick',[0 25 50])
set(gca,'position',[0.1954    0.1700    0.7096    0.27])
grid on
% subplot(2,1,2)
% semilogy([1 2 4 6 8 10 15 20 25 30 35 40 45 50 55],Nd_1(2:end),'k','linewidth',3); hold on
% semilogy([1 2 4 6 8 10 15 20 25 30 35 40 45 50 55],Nd_10(2:end),'r','linewidth',3); hold on
% set(gca,'ylim',[1,1e12],'fontsize',14,'ytick',10.^(0:6:12),'xlim',[0 55],'xtick',[0 25 50])
% set(gca,'position',[0.1954    0.1700    0.7096    0.22])
% grid on
subplot(2,1,1)
semilogy([1 2 4 6 8 10 15 20 25 30 35 40 45 50 55],Nm_1(2:end),'k','linewidth',3); hold on
semilogy([1 2 4 6 8 10 15 20 25 30 35 40 45 50 55],Nm_10(2:end),'r','linewidth',3); hold on
set(gca,'ylim',[1e2,1e8],'fontsize',14,'ytick',10.^(2:2:8),'xlim',[0 55],'xtick',[0 25 50],'xticklabel',[])
set(gca,'position',[0.1954    0.5838    0.7096    0.27])
grid on
% 

%%
% figure(14)
% Nsum = 1e9;
% B = linspace(sig2true*0.99,sig2true*1.01,1000);
% Xi2_range = (Nsum-1)*B/sig2true;
% Q = chi2pdf(Xi2_range,Nsum-1);
% plot(B,Q/sum(Q),'c','linewidth',2); hold on
% 
% s1111 = P*(([0:N]-mn).^4)';
% mu4 = s1111-(Nsum-3)/(Nsum-1)*sig2true^2;
% Q = normpdf(B,sig2true,sqrt(1/Nsum*mu4));
% plot(B,Q/sum(Q),'m','linewidth',2)
% title([num2str(Nsum),' cells']);
