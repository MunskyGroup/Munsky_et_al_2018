function Make_Panels_For_Figure_2(GENE)
%% This will make figure 3.
SPACE = 'NonSpatial';  %'Spatial' or 'NonSpatial'
% GENE = 'STL1';      %'STL1' or 'CTT1'
DIST = {'Means','Moments_Corrected','Distributions','Moments_4'};
path(path,'Solvers/moments')
load(['Fit_Results/',GENE,'_0p2and0p4/NonSpatial/Moments_4/4State_k21only_noloop/Analyses.mat']);

%% Plot Means
f1 = figure(21);clf;
f121 = figure(121);clf;
set(f1,'position',[1000        1194         607         144]);
for k=1:4
    try close(1)
    catch
    end
    uiopen(['Figures/',GENE,'_',SPACE,'_Fit_',DIST{k},'_Plot_Means.fig'],1)
    figure(1)
    clear TF
    for j=1:2
        subplot(1,2,j)
        TF(j).ax = gca;
    end
    
    if k<4
        figure(f1);subplot(1,3,k)
    else
        figure(f121);
    end
    
    for j=1:2
        copyobj(allchild(TF(j).ax),gca);
    end
    set(gca,'fontsize',15,'xticklabel',[],'ylim',[0 80])
    
    pla = get(gca,'Children');
    set(pla(1:3),'markersize',10,'marker','^')
    set(pla(5:7),'markersize',10,'marker','o')
    
    if k==1||k==4
    else
        set(gca,'yticklabel',[])
    end
    
    if k<4
        P = get(gca,'position');
        set(gca,'position',[P(1)-0.05,P(2)-0.05,0.25 0.8])
    else
        FUNNAME = ['conv_mom_1_4'];
        funhandle = str2func(FUNNAME);
        load(['Fit_Results/',GENE,'_0p2and0p4/NonSpatial/Moments_4/4State_k21only_noloop/Analyses.mat']);
        inds = find(max(Hog.MOD.Data(1).Bred(2:end,1:4),[],2)==0);
        Hog.MOD = struct(Hog.MOD);
        for iss = 1:2
            for i_time = 1:16
                n = Hog.MOD.Data(iss).Num_cells(i_time);
                xm_uncen = Hog.Moments_4(iss).moments.Trajectories(i_time,inds)';
                for i=length(xm_uncen):-1:1; v{i} = xm_uncen(i); end
                xm = funhandle(v{:});
                VL(i_time) = xm(2)/n;
                mL(i_time) = Hog.Moments_4(iss).moments.Trajectories(i_time,inds(1))';
                %                 EM(i_time) = (xm(4) - (n-3)/(n-1)*(xm(2))^2);
            end
            hold on;
            if iss==1
                H = errorbar(Hog.MOD.Data(iss).Times/60,mL,sqrt(VL),'r','linewidth',2);
            elseif iss==2
                H = errorbar(Hog.MOD.Data(iss).Times/60+0.5,mL,sqrt(VL),'b','linewidth',2);
            end
            set(gca,'ylim',[0,400],'xtick',[0:10:60],'xticklabel',[0:10:60])
        end

    end
    box on
end


%% Plot Standard Deviations
f2 = figure(22);clf;
f122 = figure(122);clf;
set(f2,'position',[1000        994         607         144]);
for k=1:4
    try close(1)
    catch
    end
    uiopen(['Figures/',GENE,'_',SPACE,'_Fit_',DIST{k},'_Plot_Vars.fig'],1)
    figure(1)
    for j=1:2
        subplot(1,2,j)
        TF(j).ax = gca;
    end
    
    if k<4
        figure(f2);subplot(1,3,k)
    else
        figure(f122);
    end
    for j=1:2
        copyobj(allchild(TF(j).ax),gca);
    end
    set(gca,'fontsize',15,'xticklabel',[],'ylim',[0 60])
    
    pla = get(gca,'Children');
    set(pla(1:3),'markersize',10,'marker','^')
    set(pla(5:7),'markersize',10,'marker','o')
    
    if k==1||k==4
    else
        set(gca,'yticklabel',[])
    end
    
    if k<4
    P = get(gca,'position');
    set(gca,'position',[P(1)-0.05,P(2)-0.05,0.25 0.8])
    else
        for kk=[1,2,3,5,6,7]
            pla(kk).YData=pla(kk).YData.^2;
        end
        
        FUNNAME = ['conv_mom_1_4'];
        funhandle = str2func(FUNNAME);
        load(['Fit_Results/',GENE,'_0p2and0p4/NonSpatial/Moments_4/4State_k21only_noloop/Analyses.mat']);
        inds = find(max(Hog.MOD.Data(1).Bred(2:end,1:4),[],2)==0);
        Hog.MOD = struct(Hog.MOD);
        for iss = 1:2
            for i_time = 1:16
                n = Hog.MOD.Data(iss).Num_cells(i_time);
                xm_uncen = Hog.Moments_4(iss).moments.Trajectories(i_time,inds)';
                for i=length(xm_uncen):-1:1; v{i} = xm_uncen(i); end
                xm = funhandle(v{:});
                VL(i_time) = xm(2);
                EM(i_time) = (xm(4) - (n-3)/(n-1)*(xm(2))^2)/n;
            end
            hold on;
            UEB = EM.^(1/2);
            LEB = min(VL-1e-6,UEB);
            if iss==1
                H = errorbar(Hog.MOD.Data(iss).Times/60,(VL),LEB,UEB,'r','linewidth',2)
            elseif iss==2
                H = errorbar(Hog.MOD.Data(iss).Times/60+0.5,(VL),LEB,UEB,'b','linewidth',2)
            end
            set(gca,'ylim',[0,1500],'xtick',[0:10:60],'xticklabel',[0:10:60])
        end
        
    end
    box on
    close(1)
end


%% Plot ON fractions
f3 = figure(23);clf;
set(f3,'position',[1000        794         607         144]);
for k=1:3
    try close(1)
    catch
    end
    uiopen(['Figures/',GENE,'_',SPACE,'_Fit_',DIST{k},'_Plot_ON_Total.fig'],1)
    figure(1)
    for j=1:2
        subplot(1,2,j)
        TF(j).ax = gca;
    end
    
    figure(f3);subplot(1,3,k)
    for j=1:2
        copyobj(allchild(TF(j).ax),gca);
    end
    set(gca,'fontsize',15,'ylim',[0 1.05],'ytick',[0,0.2,0.4,0.6,0.8,1],'yticklabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})
    xlabel('Time (min)')
    
    pla = get(gca,'Children');
    set(pla(1:3),'markersize',10,'marker','^')
    set(pla(5:7),'markersize',10,'marker','o')
    
    if k==1
    else
        set(gca,'yticklabel',[])
    end
    P = get(gca,'position');
    set(gca,'position',[P(1)-0.05,P(2)-0.05,0.25 0.8])
    box on
    close(1)
end

%% Plot Distributions
f4 = figure(24);clf;
set(f4,'position',[1000        594         607         144]);
time_point = 8;
iSALT=1;
for k=1:3
    if iSALT==2
        uiopen(['Figures/',GENE,'_',SPACE,'_Fit_',DIST{k},'_Plot_Marg_Dists_0p4.fig'],1);
    elseif iSALT==1
        uiopen(['Figures/',GENE,'_',SPACE,'_Fit_',DIST{k},'_Plot_Marg_Dists_0p2.fig'],1);
    end
    f1= figure(1);
    
    figure(f4);subplot(1,3,k)
    copyobj(allchild(f1.Children(time_point)),gca);
    
    set(gca,'fontsize',15,'ylim',[0 0.05],'xlim',[0 160],'ytick',[0.00,0.02,0.04],'yticklabel',{'0.00','0.02','0.04'})
    xlabel('Number of mRNA')
    pla =get(gca,'children');
    set(pla(1:3),'Visible','off')
    set(pla(4:6),'linewidth',3)
    set(pla(7),'linewidth',5)
    
    
    if k==1
%         ylabel('Probability')
    else
        set(gca,'yticklabel',[])
    end
    P = get(gca,'position');
    set(gca,'position',[P(1)-0.05,P(2)+0.05,0.25 0.8])
    box on
    close(1)
end

