function Solve_GUI(GUI)
%% Moments
if GUI.PlotVariancesButton.Value==1||GUI.PlotMeansButton.Value==1
    if isempty(GUI.Moments)
        disp('Solving moments');
        tic
        GUI.MOD.Salt = [0.2 0.4];
        GUI.MOD.Type='Moments';
        GUI.MOD.Fit_Type='Moments';
        GUI.Moments = Hog_Model.solve_hog(GUI.MOD);
        disp(['Done in ',num2str(toc),'s.']);
    end
    
    Hog_Model.plot_moments_mod(GUI.MOD,GUI.Moments,...
        GUI.FigEditField.Value,GUI.FigEditField_2.Value);
    for TMP = {'_A','_B','_C'}
        GUI.MOD.Reps = TMP{1};
        Hog_Model.plot_moments_dat(GUI.MOD,GUI.MOD.Data,...
            GUI.FigEditField.Value,GUI.FigEditField_2.Value);
    end
    f1=figure(GUI.FigEditField.Value);
    f2=figure(GUI.FigEditField_2.Value);
    Adjust_Plots.Mean(f1,GUI.MOD);
    Adjust_Plots.Vars(f2,GUI.MOD);
end
%% Distributions
if GUI.PlotDistributionsButton.Value==1
    if isempty(GUI.Distributions)
        disp('Solving distributions');
        tic
        GUI.MOD.Salt = [0.2 0.4];
        GUI.MOD.Type='Distributions';
        GUI.Distributions = Hog_Model.solve_hog(GUI.MOD);
        disp(['Done in ',num2str(toc),'s.']);
    end
    
    Hog_Model.plot_distributions_mod(GUI.MOD,GUI.Distributions,...
        GUI.FigEditField_3.Value,GUI.FigEditField_3.Value+2);
    for TMP = {'_A','_B','_C'}
        GUI.MOD.Reps = TMP{1};
        Hog_Model.plot_distributions_dat(GUI.MOD,GUI.MOD.Data,...
            GUI.FigEditField_3.Value,GUI.FigEditField_3.Value+2);
    end
    Adjust_Plots.Marginal(GUI.FigEditField_3.Value,GUI.MOD,0.2)
    Adjust_Plots.Marginal(GUI.FigEditField_3.Value+1,GUI.MOD,0.4)
    Adjust_Plots.On(GUI.FigEditField_3.Value+2,'NonSpatial')
    
    switch GUI.SpatialPlotButtonGroup.SelectedObject.Text
        case 'Spatial'
            GUI.MOD.Type='Distributions';
            Hog_Model.plot_spat_distributions_mod(GUI.MOD,GUI.Distributions,...
                1000+GUI.FigEditField_3.Value,1000+GUI.FigEditField_3.Value+4);
            for TMP = {'_A','_B','_C'}
                GUI.MOD.Reps = TMP{1};
                Hog_Model.plot_spat_distributions_dat(GUI.MOD,GUI.MOD.Data,...
                    1000+GUI.FigEditField_3.Value+2,1000+GUI.FigEditField_3.Value+4);
            end
            Adjust_Plots.On(1000+GUI.FigEditField_3.Value+4,'Spatial')
            
            GUI.MOD.Reps = '_AC';
            Hog_Model.plot_spat_distributions_dat(GUI.MOD,GUI.MOD.Data,...
                1000+GUI.FigEditField_3.Value+2,500); close(500);
            Adjust_Plots.Joint(1000+GUI.FigEditField_3.Value,GUI.MOD,0.2,'Model','Cyt','Nuc')
            Adjust_Plots.Joint(1000+GUI.FigEditField_3.Value+1,GUI.MOD,0.4,'Model','Cyt','Nuc')
            Adjust_Plots.Joint(1000+GUI.FigEditField_3.Value+2,GUI.MOD,0.2,'Data','Cyt','Nuc')
            Adjust_Plots.Joint(1000+GUI.FigEditField_3.Value+3,GUI.MOD,0.4,'Data','Cyt','Nuc')
        case 'NonSpatial'
            for f=1237:1238
                for v = 1:5%1:length(obj.tt)
                    COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4];
                    W = get(f,'Children');
                    C = get(W(end+1-v),'Children');
                    mn = C(1).YData;
                    mx = C(1).YData;
                    for i=1:3
                        set(C(i),'Color',COL(i,:))
                        TMP = length(C(i).YData);
                        if length(mn)<TMP
                            mn(length(mn)+1:TMP) = 0;
                            mx(length(mx)+1:TMP) = C(i).YData(length(mx)+1:TMP);
                        end
                        mn = min(mn(1:TMP),C(i).YData(1:TMP));
                        mx = max(mx(1:TMP),C(i).YData(1:TMP));
                    end
                    
                    Q = figure(f+2);
                    set(Q,'position',[514   754   578    91]);
                    switch GUI.AnalysisTypeButtonGroup.SelectedObject.Text
                        case 'FSP'
                            cl = 'k';
                        case 'Extended Moments'
                            cl = 'm';
                        case '2nd Moments'
                            cl = 'b';
                        case 'Means'
                            cl = 'r';
                    end
                       
                    if length(Q.Children)<v
                        mn=max(1e-5,mn);
                        mx=max(1e-5,mx);
                        subplot(1,5,v)
                        set(gca,'position',[0.05+0.13*(v-1) 0.1100 0.1237 0.8150],'yticklabel',[],'xticklabel',[])
                        fill([0:length(mn)-1,length(mn)-1:-1:0,0],[mn(1:end),mx(end:-1:1),mn(1)],[0.6 0.6 0.6]); hold on
                        set(gca,'xlim',[-5 100],'ylim',[0 0.04],'yticklabel',[],'xticklabel',[],'yscale','lin');
                        pp=plot([0:length(C(4).YData)-1],C(4).YData,cl,'linewidth',3);
                        pp.Color(4) = 0.5;
                    else
                        pp = plot(Q.Children(end+1-v),[0:length(C(4).YData)-1],C(4).YData,cl,'linewidth',3);
                        pp.Color(4) = 0.5;
                    end
                end
            end
            TPR = 0.7;
            %%
            Q1 = figure(203);
            for i=1:2
                C = Q1.Children(2:-1:1);
                x = C(i).Children(1).XData;
                mn = C(i).Children(1).YData;
                mx = C(i).Children(1).YData;
                for i2=1:3
                    mn = min(mn,C(i).Children(i2).YData);
                    mx = max(mx,C(i).Children(i2).YData);
                end
                Q2 = figure(205);
                set(Q2,'Name','On vs. time')
                set(Q2,'Position',[1000         626         205         329]);
                subplot(2,1,i)
                if length(Q2.Children(1).Children)<=1
                    fill([x(1:end),x(end:-1:1),x(1)],[mn(1:end),mx(end:-1:1),mn(1)],[0.6 0.6 0.6]); hold on
                    set(gca,'xticklabel',[],'yticklabel',[])
                end
                pp = plot(x,C(i).Children(4).YData,cl,'linewidth',3); hold on;
                pp.Color(4) = TPR;
            end
            %%
            if GUI.NonSpatialButton_2.Value
                Q1 = figure(1);
                for i=1:2
                    C = Q1.Children(2:-1:1);
                    x = C(i).Children(1).XData;
                    mn = C(i).Children(1).YData;
                    mx = C(i).Children(1).YData;
                    for i2=1:3
                        mn = min(mn,C(i).Children(i2).YData);
                        mx = max(mx,C(i).Children(i2).YData);
                    end
                    Q2 = figure(206);
                    set(Q2,'Name','Mean vs. time')
                    set(Q2,'Position',[1000         626         205         329]);
                    subplot(2,1,i)
                    if length(Q2.Children(1).Children)<=1
                        fill([x(1:end),x(end:-1:1),x(1)],[mn(1:end),mx(end:-1:1),mn(1)],[0.6 0.6 0.6]); hold on
                        set(gca,'xticklabel',[],'yticklabel',[])
                    end
                    pp = plot(x,C(i).Children(4).YData,cl,'linewidth',3); hold on;
                    pp.Color(4) = TPR;
                end
                %%
                Q1 = figure(101);
                for i=1:2
                    C = Q1.Children(2:-1:1);
                    x = C(i).Children(1).XData;
                    mn = C(i).Children(1).YData;
                    mx = C(i).Children(1).YData;
                    for i2=1:3
                        mn = min(mn,C(i).Children(i2).YData);
                        mx = max(mx,C(i).Children(i2).YData);
                    end
                    Q2 = figure(207);
                    set(Q2,'Name','Var vs. time')
                    set(Q2,'Position',[1000         626         205         329]);
                    subplot(2,1,i)
                    if length(Q2.Children(1).Children)<=1
                        fill([x(1:end),x(end:-1:1),x(1)],[mn(1:end),mx(end:-1:1),mn(1)],[0.6 0.6 0.6]); hold on
                        set(gca,'xticklabel',[],'yticklabel',[])
                    end
                    pp = plot(x,C(i).Children(4).YData,cl,'linewidth',3); hold on;
                    pp.Color(4) = TPR;
                end
            end

            if strcmp(cl,'k')
                %%
                for IL =  205:207
                    Q2 = figure(IL);
                    if ~isempty(Q2.Children)
                        set(Q2.Children(2),'Position',[0.1300    0.51    0.77    0.34]);
                        set(Q2.Children(1),'Position',[0.1300    0.11    0.77    0.34]);
                        switch IL
                            case 205
                                set(Q2,'Position',[1000         627         154         328]);
                                set(Q2.Children(1),'ylim',[0 1])
                                set(Q2.Children(2),'ylim',[0 1])
                                for k=1:length(Q2.Children(1).Children)-1
                                    Q2.Children(1).Children(k).Color(4) = 0.7;
                                    Q2.Children(2).Children(k).Color(4) = 0.7;
                                end
                            case 206
                                set(Q2,'Position',[1000         627         154         328]);
                                set(Q2.Children(1),'ylim',[0 80])
                                set(Q2.Children(2),'ylim',[0 80])
                                for k=1:length(Q2.Children(1).Children)-1
                                    Q2.Children(1).Children(k).Color(4) = 0.7;
                                    Q2.Children(2).Children(k).Color(4) = 0.7;
                                end
                            case 207
                                set(Q2,'Position',[1000         627         154         328]);
                                set(Q2.Children(1),'ylim',[0 80])
                                set(Q2.Children(2),'ylim',[0 80])
                                for k=1:length(Q2.Children(1).Children)-1
                                    Q2.Children(1).Children(k).Color(4) = 0.7;
                                    Q2.Children(2).Children(k).Color(4) = 0.7;
                                end
                        end
                    end
                end
            end
    end
end
%% Transcription sites.
if GUI.NascentmRNAButton.Value==1
    TS_Plots_for_GUI
    
    figure(313)
    set(gcf,'position',[601   698   180   110])
    set(gcf,'Name','TS Fraction versus time')
    q = get(gca,'children');
    set(q([1,3]),'MarkerSize',12,'Marker','x');
    set(q(2),'color','k');
    set(gca,'xlim',[0,55],'xtick',[0 25 50])
    
    figure(323)
    set(gcf,'position',[535   455  180   110])
    set(gcf,'Name','TS Fraction versus time')
    q = get(gca,'children');
    set(q([1,2,4]),'MarkerSize',12,'Marker','x');
    set(q(3),'color','k');
    set(gca,'xlim',[0,55],'xtick',[0 25 50])
    
    figure(314)
    set(gcf,'Name','Nascent mRNA Distributions versus time (0.2M NaCl)')
    figure(324)
    set(gcf,'Name','Nascent mRNA Distributions versus time (0.4M NaCl)')
    
    figure(304)
    xlabel('Number mRNA')
    ylabel('Probability')
    set(gcf,'Name','Time average of Nascent mRNA distribution (0.2M NaCl)')

    figure(305)
    xlabel('Number mRNA')
    ylabel('Probability')
    set(gcf,'Name','Time average of Nascent mRNA distribution (0.4M NaCl)')
    
    Q2 = figure(1239);
    if ~isempty(Q2.Children)
        set(gcf,'Name','Nascent mRNA Distributions vs time (0.2M NaCl)')
    else
        close(Q2);
    end
    Q2 = figure(1240);
    if ~isempty(Q2.Children)
        set(gcf,'Name','Nascent mRNA Distributions vs time (0.4M NaCl)')
    else
        close(Q2);
    end
end
%%
Q2 = figure(1237);
if ~isempty(Q2.Children)
    set(gcf,'Name','Mature mRNA Distributions vs Time (0.2M NaCl')
else
    close(Q2);
end
Q2 = figure(1238);
if ~isempty(Q2.Children)
    set(gcf,'Name','Mature mRNA Distributions vs Time (0.4M NaCl')
else
    close(Q2);
end

for j=339:340
Q = figure(j);
for i=1:5
    try
        switch j
            case 339
                set(gcf,'Name','Nascent mRNA Distributions vs time (0.2M NaCl)')
            case 340
                set(gcf,'Name','Nascent mRNA Distributions vs time (0.4M NaCl)')
        end
        set(Q.Children(i),'position',[0.0300+0.13*(6-i)    0.1650    0.1237    0.7600],...
            'xticklabel',[]);
    catch
    end
end
end
