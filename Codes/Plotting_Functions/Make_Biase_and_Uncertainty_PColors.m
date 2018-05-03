% PAR_NAMES = {'k_{12}','k_{21}','k_{23}','k_{32}','k_{34}','k_{43}','k_{21}','k_{r1}','k_{r2}','k_{r3}','k_{r4}','\gamma','t_{offset}'};
clear all
PAR_NAMES_Cell = {'k_{12}','k_{23}','k_{32}','k_{34}','k_{43}','k_{r2}','k_{r3}','k_{r4}','\gamma'};

% figure(1);
% set(1,'position',[ 242        1116         781         178]);
% pcolor([Hog.MC_Results.Bias_Indv1,zeros(3,1);zeros(1,length(PAR_NAMES)+1)])
% set(gca,'fontsize',16,'xtick',[1.5:1:length(PAR_NAMES)+.5],'xticklabel',PAR_NAMES)
% set(gca,'ytick',[1.5:1:3.5],'yticklabel',{'\mu','\mu,\sigma','FSP'})
% set(gca,'XAxisLocation','top')
% cb = colorbar;
% set(cb.Label,'String','log_{10}(bias)')
%
%%
close all
for gg = 1:2
    switch gg
        case 2
            gene = 'STL1';
        case 1
            gene = 'CTT1';
    end
    load(['TMP_MC_',gene,'.mat']);
    %     for i=[1 2 3 4 5 6]
    for i=1:10
        switch i
            case 1
                ty = 'MeansNonSpatial';
                ttl = [gene,' Means NonSpatial'];
            case 2
                ty = 'Moments_CorrectedNonSpatial';
                ttl = [gene,' Moments NonSpatial'];
            case 3
                ty = 'DistributionsNonSpatial';
                ttl = [gene,' Distributions NonSpatial'];
            case 4
                ty = 'Moments_4NonSpatial';
                ttl = [gene,' First 4 Moments NonSpatial'];
            case 5
                ty = 'GaussianNonSpatial';
                ttl = [gene,' Gaussian'];
            case 6
                ty = 'MeansSpatial';
                ttl = [gene,' Means Spatial'];
            case 7
                ty = 'Moments_CorrectedSpatial';
                ttl = [gene,' Moments Spatial'];
            case 8
                ty = 'DistributionsSpatial';
                ttl = [gene,' Distributions Spatial'];
            case 9
                ty = 'Moments_4Spatial';
                ttl = [gene,' First 4 Moments Spatial'];
            case 10
                ty = 'GaussianSpatial';
                ttl = [gene,' Gaussian Spatial'];
        end
        figure(1+i+(gg-1)*10); clf;
        set(1+i+(gg-1)*8,'position',[8+800*(mod(i-1,4))   80+600*floor((i-1)/4)   781   502]);
        
        subplot(2,1,1)
        rng_low = 0;
        rng_high = 5;
        TMP = max((Hog.MC_Results.Bias_Indv1+Hog.MC_Results.Bias_Indv2)/2,rng_low);
        TMP = min(TMP,rng_high);
        pcolor([TMP(i,:),rng_low*ones(1,1);rng_high*ones(1,length(PAR_NAMES_Cell)+1)])
        
        
        set(gca,'fontsize',16,'xtick',[1.5:1:length(PAR_NAMES_Cell)+.5])
        if i==1||i==5
            set(gca,'xticklabel',PAR_NAMES_Cell)
            set(gca,'XAxisLocation','top')
        else
            set(gca,'xticklabel',[]);
        end
        
        set(gca,'ytick',[])
        cb1 = colorbar;
        set(gca,'position',[0.13      0.7      0.71098      0.15])
        %         set(cb1.Label,'String','log_{10}(bias)')
        title(ttl)
        for ii=1:length(PAR_NAMES_Cell)
            if Hog.MC_Results.Bias_Sign(i,ii)==1
                tx = '+';
            else
                tx = '-';
            end
            t = text(ii+0.5,1.5,tx);
            set(t,'fontsize',30,'color',[0.9 0.9 0.9])
        end
        
        subplot(2,1,2)
        A = getfield(Hog.MC_Results,ty);
        rng_low = -10;
        rng_high = 5;
        TMP = [A.COV_log(end:-1:1,:),rng_low*ones(length(PAR_NAMES_Cell),1);rng_high*ones(1,length(PAR_NAMES_Cell)+1)];
        TMP = max(TMP,rng_low);
        TMP = min(TMP,rng_high);
        pcolor(TMP)
        
        set(gca,'fontsize',16,'xtick',[1.5:1:length(PAR_NAMES_Cell)+.5],'ytick',[1.5:1:length(PAR_NAMES_Cell)+.5])
        
        if i==1||i==2||i==3||i==4||i==5
            set(gca,'yticklabel',PAR_NAMES_Cell(end:-1:1))
        else
            set(gca,'yticklabel',[]);
        end
        if i==3||i==8
            set(gca,'xticklabel',PAR_NAMES_Cell)
            set(gca,'XAxisLocation','bottom');
        else
            set(gca,'xticklabel',[]);
        end
        
        cb2 = colorbar;
        %         set(cb2.Label,'String','log_{10}(|\sigma_{ij}|)')
        %     colormap('summer')
        
        for ii=1:length(PAR_NAMES_Cell)
            for jj=ii:length(PAR_NAMES_Cell)
                if A.COV_sign(length(PAR_NAMES_Cell)-ii+1,length(PAR_NAMES_Cell)-jj+1)==1
                    tx = '+';
                else
                    tx = '-';
                end
                t = text(-ii+1.5+length(PAR_NAMES_Cell),jj+0.5,tx);
                set(t,'fontsize',30,'color',[0.9 0.9 0.9])
                t = text(-jj+1.5+length(PAR_NAMES_Cell),ii+0.5,tx);
                set(t,'fontsize',30,'color',[0.9 0.9 0.9])
                
            end
        end
        set(gca,'position',[0.13      0.1      0.71098      0.55])
    end
end


