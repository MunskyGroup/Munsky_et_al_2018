function Uncrt = Update_MCMC_Plots
MOD = Hog_Model;
MOD.Recompile_MCMC = 0;
fnum=20;
% Par_set_ars = {'transit_pars','All_Pars','kr3_and_gamma'};
  Par_set_ars = {'All_Pars','Spat_Pars','Red_Pars'};
%  Par_set_ars = {'All_Pars'};
% Par_set_ars = {'Red_Pars'};
comp_gene = 'STL1';

for parsets = Par_set_ars
    switch parsets{1}
        case 'All_Pars'
            PARINDS_SP = [1:11,14,15]; PARINDS_NSP =[1:13]; inset = [10 12]; scaling=1;
            GENES_AR = {'CTT1','STL1'};
%             PAR_NAMES = 'Gene,Stats,Spatial,k12a,k21a,k23a,k32a,k34a,k43a,k21b,kr1,kr2,kr3,kr4,gamma,t_offset';
            PAR_NAMES_b = 'Gene,&,Stats,&,Spatial,&,k12a,&,std,&,k21a,&,std,&,k23a,&,std,&,k32a,&,std,&,k34a,&,std,&,k43a,&,std,&,k21b,&,std,&,kr1,&,std,&,kr2,&,std,&,kr3,&,std,&,kr4,&,std,&,gamma,&,std,&,t_offset,&,std';
            MOD.MCMC_Reps=0;
        case 'Red_Pars'
            PARINDS_SP = [1,3:6,9:11,14]; PARINDS_NSP =[1,3:6,9:11,12]; inset = [7 9]; scaling=1;
            GENES_AR = {'CTT1','STL1'};
%             PAR_NAMES = 'Gene,Stats,Spatial,k12a,k23a,k32a,k34a,k43a,kr2,kr3,kr4,gamma';
            PAR_NAMES_b = 'Gene,&,Stats,&,Spatial,&,k12a,&,std,&,k23a,&,std,&,k32a,&,std,&,k34a,&,std,&,k43a,&,std,&,kr2,&,std,&,kr3,&,std,&,kr4,&,std,&,gamma,&,std';
            MOD.MCMC_Reps=1;
        case 'Spat_Pars'
            PARINDS_SP = [12,14,15]; inset = [2 3]; scaling=1;
            GENES_AR = {'CTT1','STL1'};
%             PAR_NAMES = 'Gene,Stats,Spatial,gamma_n,gamma_c,transport';
            PAR_NAMES_b = 'Gene,&,Stats,&,Spatial,&,gamma_n,&,std,&,gamma_c,&,std,&,transport,&,std';
            MOD.MCMC_Reps=1;
        case 'kr3_and_gamma'
            PARINDS_SP = [10,14]; PARINDS_NSP =[10,12]; inset = []; scaling=0;
            GENES_AR = {'CTT1','STL1'};
%             PAR_NAMES = 'Gene,Stats,Spatial,kr3,gamma';
            PAR_NAMES_b = 'Gene,&,Stats,&,Spatial,&,kr3,&,std,&,gamma,&,std';
            MOD.MCMC_Reps=1;
        case 'transit_pars'
            PARINDS_SP = [1:7]; PARINDS_NSP =[1:7];  inset = [3 5]; scaling=1;
%             GENES_AR = {'STL1','CTT1','STL1_and_CTT1'};
            GENES_AR = {'CTT1','STL1'};
%             PAR_NAMES = 'Gene,Stats,Spatial,k12a,k21a,k23a,k32a,k34a,k43a,k21b';
            PAR_NAMES_b = 'Gene,&,Stats,&,Spatial,&,k12a,&,std,&,k21a,&,std,&,k23a,&,std,&,k32a,&,std,&,k34a,&,std,&,k43a,&,std,&,k21b,&,std';
            MOD.MCMC_Reps=1;
    end
    fid = fopen(['FOUND_PARS_',parsets{1},'.csv'],'w');
    fprintf(fid,PAR_NAMES_b);
    Hog.MC_Results.Bias = [];
    Hog.MC_Results.Bias_Sign = [];
    Hog.MC_Results.Bias_Indv1 = [];
    Hog.MC_Results.Bias_Indv2 = [];
    Hog.MC_Results.Uncert = [];

    for GENES = GENES_AR
        Hog.MOD = MOD;
        Hog.MOD.Gene = GENES{1};
        switch GENES{1}
            case {'STL1','CTT1'}
                Hog.MOD.Model = '4State_k21only_noloop';
                TYPE_AR = {'Means','Moments_Corrected','Distributions','Moments_4','Gaussian'};
%                 TYPE_AR = {'Means','Moments_Corrected','Distributions'};
                DIR = ['Fit_Results/',Hog.MOD.Gene,'_0p2and0p4/Spatial/Distributions/',Hog.MOD.Model]
                Hog.MOD.Spatial = 'Spatial';
                Hog.MOD.Fit_Type = 'Distributions';
                Hog.MOD.MCMC_File = [DIR,'/MCMC_Processed.mat'];
                BEST = log10(abs(Hog.MOD.MCMC_Chain_Obj.Best(PARINDS_SP)));
                fnum=20;
                Hog.MC_Results.Bias = [];
                Hog.MC_Results.Bias_Sign = [];
                Hog.MC_Results.Uncert = [];
                Hog.MC_Results.Bias_Indv1 = [];
                Hog.MC_Results.Bias_Indv2 = [];
            case 'STL1_and_CTT1'
                Hog.MOD.Model = '4State_2gene_dep_k21only_noloop';
                TYPE_AR = {'Means','Moments_Corrected','Moments_4'};
                DIR = ['Fit_Results/',comp_gene,'_0p2and0p4/Spatial/Distributions/4State_k21only_noloop'];
                Hog.MOD.Spatial = 'Spatial';
                Hog.MOD.Fit_Type = 'Distributions';
                Hog.MOD.MCMC_File = [DIR,'/MCMC_Processed.mat'];
                BEST = log10(abs(Hog.MOD.MCMC_Chain_Obj.Best(PARINDS_SP)));
                Hog.MC_Results.Bias_Sign = [];
                Hog.MC_Results.Bias = [];
                Hog.MC_Results.Uncert = [];
                Hog.MC_Results.Bias_Indv1 = [];
                Hog.MC_Results.Bias_Indv2 = [];
        end
        
        if ~strcmp(parsets{1},'Spat_Pars')
            SP_Ar = {'NonSpatial','Spatial'};
        else
            SP_Ar = {'Spatial'};
        end 
        for SP=SP_Ar
            Hog.MOD.Spatial=SP{1};
            close all
            switch SP{1}
                case 'Spatial'
                    PARINDS=PARINDS_SP;
                    g=2;
                case 'NonSpatial'
                    PARINDS=PARINDS_NSP;
                    g=1;
            end
            
            if strcmp(GENES{1},'STL1_and_CTT1')
                fnum = open(['Figures/MCMC_',comp_gene,'_',SP{1},'_',parsets{1},'.fig']);
                figure(fnum); hold on
            end
            
            for FT = TYPE_AR
                Hog.MOD.Fit_Type = FT{1};
                switch GENES{1}
                    case {'STL1','CTT1'}
                        DIR = ['Fit_Results/',Hog.MOD.Gene,'_0p2and0p4/',Hog.MOD.Spatial,'/',Hog.MOD.Fit_Type,'/',Hog.MOD.Model];
                    case 'STL1_and_CTT1'
                        DIR = ['Fit_Results/',Hog.MOD.Gene,'/',Hog.MOD.Spatial,'/',Hog.MOD.Fit_Type,'/',Hog.MOD.Model];
                end
                Hog.MOD.MCMC_File = [DIR,'/MCMC_Processed.mat'];
                MCMC_Results = Hog_Model.plot_MCMC_results(Hog.MOD,Hog.MOD.MCMC_Chain_Obj,PARINDS,fnum);
                
                fprintf(fid,['\n ',GENES{1},',&,',FT{1},',&,',SP{1},',&,']);
                for k=1:length(MCMC_Results.M);
                    fprintf(fid,num2str(MCMC_Results.M(k)));fprintf(fid,',&,');
                    fprintf(fid,num2str(MCMC_Results.STD(k)));fprintf(fid,',&,');
                
                end
                
                Hog.MC_Results.Uncert = [Hog.MC_Results.Uncert;...
                    MCMC_Results.trace,MCMC_Results.trace1,MCMC_Results.trace2];
                
                Hog.MC_Results.Bias = [Hog.MC_Results.Bias;...
                            norm(abs(log10(abs(MCMC_Results.M./10.^BEST)))),...
                            norm(abs(log10(abs(MCMC_Results.M1./10.^BEST)))),...
                            norm(abs(log10(abs(MCMC_Results.M2./10.^BEST))))];

                Hog.MC_Results.Bias_Indv1 = [Hog.MC_Results.Bias_Indv1;...
                            abs(log10(abs(MCMC_Results.M1./10.^BEST)))];
                Hog.MC_Results.Bias_Indv2 = [Hog.MC_Results.Bias_Indv2;...
                            abs(log10(abs(MCMC_Results.M2./10.^BEST)))];
                        
                Hog.MC_Results.Bias_Sign = [Hog.MC_Results.Bias_Sign;...
                            sign(MCMC_Results.M-10.^BEST)];
                
                Hog.MC_Results = setfield(Hog.MC_Results,[FT{1},SP{1}],MCMC_Results); 
                
                
                if strcmp(parsets{1},'All_Pars')
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).STD = MCMC_Results.STD;
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).M = MCMC_Results.M;
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).UB99 = MCMC_Results.UB99;
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).LB99 = MCMC_Results.LB99;
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).UB90 = MCMC_Results.UB90;
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).LB90 = MCMC_Results.LB90;
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).UB50 = MCMC_Results.UB50;
                    Uncrt.(GENES{1}).(FT{1}).(SP{1}).LB50 = MCMC_Results.LB50;
                end
            end
            
            
            FIG_Sub = ['Figures/MCMC_',GENES{1},'_',SP{1},'_',parsets{1},'.fig'];
            saveas(fnum,FIG_Sub,'fig');
            Adjust_MCMC_Plots(fnum,PARINDS,scaling,BEST,Hog.MOD,inset,0,MOD.MCMC_Reps)
            FIG_Sub = ['Figures/MCMC_',GENES{1},'_',SP{1},'_',parsets{1},'_Reformatted.fig'];
            saveas(fnum,FIG_Sub,'fig');
            
%             Adjust_MCMC_Plots(fnum+1,PARINDS,scaling,BEST,Hog.MOD,inset,1,1)
%             FIG_Sub = ['Figures/MCMC_',GENES{1},'_',SP{1},'_',parsets{1},'_Reformatted_Centered.fig'];
%             saveas(fnum+1,FIG_Sub,'fig');
           
        end
        if strcmp(parsets{1},'Red_Pars')
            save(['TMP_MC_',GENES{1},'.mat'])
        end

        switch GENES{1}
            case 'STL1'
                STL1_Rslts = Hog.MC_Results.Uncert;
            case 'CTT1'
                CTT1_Rslts = Hog.MC_Results.Uncert;
            case 'STL1_and_CTT1'
                Hog.MC_Results.Uncert = [STL1_Rslts;CTT1_Rslts;Hog.MC_Results.Uncert];
        end
        
        if MOD.MCMC_Reps==1
            Make_Bias_Uncertainty_Plot(Hog,3,g);
            FIG_Sub = ['Figures/MCMC_',GENES{1},'_Errors_',parsets{1},'.fig'];
            saveas(3,FIG_Sub,'fig');
            FIG_Sub2 = ['Figures/MCMC_',GENES{1},'_Errors_',parsets{1},'Indv.fig'];
            saveas(4,FIG_Sub2,'fig');
        end
    end
end

function Adjust_MCMC_Plots(f,PARINDS,scaling,BEST,obj,inset,centering,open)
if centering==1
    BEST =0*BEST;
end
COL = [0 0 0;0 0 1;1 0 0;1 0 1;0 1 1;0 1 0];
figure(f); set(f,'position',[ 347         214        1731        1063])

NP = length(PARINDS);

if exist('inset','var')&&~isempty(inset)
    subplot(NP-1,NP-1,(inset(1)-1)*(NP-1)+inset(2)-1);
    ax1=gca;
    subplot(NP-1,NP-1,(NP-1)*(NP-2)+1);
    copyobj(allchild(ax1),gca); hold on
    set(gca,'fontsize',16,'position',[0.093044     0.097497      0.35772      0.35178]);
    pla = get(gca,'children');
    if strcmp(obj.Gene,'STL1_and_CTT1')
        set(pla(1),'color',COL(4,:),'linewidth',3);
        set(pla(2),'facecolor',COL(4,:));
        set(pla(3),'color',COL(5,:),'linewidth',3);
        set(pla(4),'facecolor',COL(5,:));
        set(pla(5),'color',COL(1,:),'linewidth',3);
        set(pla(6),'facecolor',COL(1,:));
        set(pla(7),'color',COL(2,:),'linewidth',3);
        set(pla(8),'facecolor',COL(2,:));
        set(pla(9),'color',COL(3,:),'linewidth',3);
        set(pla(10),'facecolor',COL(3,:));
        set(pla(11),'color',COL(6,:),'linewidth',3);
        set(pla(12),'facecolor',COL(6,:));
        legend([pla(11),pla(9),pla(7),pla(5),pla(3),pla(1)],{'Gaussian','means','vars','dists','joint means','joint vars'})
    else
        if open==0
            set(pla(1),'color',COL(1,:),'linewidth',3);
            set(pla(2),'facecolor',COL(1,:));
            set(pla(3),'color',COL(2,:),'linewidth',3);
            set(pla(4),'facecolor',COL(2,:));
            set(pla(5),'color',COL(3,:),'linewidth',3);
            set(pla(6),'facecolor',COL(3,:));
            set(pla(7),'color',COL(4,:),'linewidth',3);
            set(pla(8),'facecolor',COL(4,:));
            set(pla(9),'color',COL(6,:),'linewidth',3);
            set(pla(10),'facecolor',COL(6,:));
            legend([pla(9),pla(7),pla(5),pla(3),pla(1)],{'means','vars','dists','vars4','Gaussian'})
        else
            set(pla(1:2),'color',COL(1,:),'linewidth',3);
            set(pla(3:4),'color',COL(2,:),'linewidth',3);
            set(pla(5:6),'color',COL(3,:),'linewidth',3);
            set(pla(7:8),'color',COL(4,:),'linewidth',3);
            set(pla(9:10),'color',COL(6,:),'linewidth',3);
%             legend([pla(5),pla(3),pla(1)],{'means','vars','dists'})
        end
    end
    xlabel(['log_{10}',obj.MCMC_Chain_Obj.par_names{PARINDS(inset(2))}]);
    ylabel(['log_{10}',obj.MCMC_Chain_Obj.par_names{PARINDS(inset(1))}]);
    for i=1:3
        YLIM = get(gca,'ylim');
        XLIM = get(gca,'xlim');
        if ~isempty(BEST)
            plot(XLIM,BEST(inset(1))*[1,1],'k--')
            plot(BEST(inset(2))*[1,1],YLIM,'k--')
        end
    end
end

for i1=1:NP-1
    for i2=i1+1:NP
        subplot(NP-1,NP-1,(i1-1)*(NP-1)+i2-1);
        set(gca,'fontsize',14);
        
        if ~isempty(BEST)
            if scaling==1
%                 YLIM = BEST(i1)+2*[-1,1];
%                 XLIM = BEST(i2)+2*[-1,1];
                XLIM(1) =min(BEST(i1)-1.1,min(get(gca,'xlim')));
                XLIM(2) =max(BEST(i1)+1.1,max(get(gca,'xlim')));
                YLIM(1) =min(BEST(i2)-1.1,min(get(gca,'ylim')));
                YLIM(2) =max(BEST(i2)+1.1,max(get(gca,'ylim')));
                set(gca,'ylim',YLIM,'xlim',XLIM,'xtick',(-100:100),'ytick',(-100:100));
                grid on
                
            end
        end
        pla = get(gca,'children');
        
    if strcmp(obj.Gene,'STL1_and_CTT1')
        set(pla(1),'color',COL(4,:),'linewidth',3);
        set(pla(2),'facecolor',COL(4,:));
        set(pla(3),'color',COL(5,:),'linewidth',3);
        set(pla(4),'facecolor',COL(5,:));
        set(pla(5),'color',COL(1,:),'linewidth',3);
        set(pla(6),'facecolor',COL(1,:));
        set(pla(7),'color',COL(2,:),'linewidth',3);
        set(pla(8),'facecolor',COL(2,:));
        set(pla(9),'color',COL(3,:),'linewidth',3);
        set(pla(10),'facecolor',COL(3,:));
        set(pla(11),'color',COL(6,:),'linewidth',3);
        set(pla(12),'facecolor',COL(6,:));
    else
        if open==0
            set(pla(1),'color',COL(1,:),'linewidth',3);
            set(pla(2),'facecolor',COL(1,:));
            set(pla(3),'color',COL(2,:),'linewidth',3);
            set(pla(4),'facecolor',COL(2,:));
            set(pla(5),'color',COL(3,:),'linewidth',3);
            set(pla(6),'facecolor',COL(3,:));
            set(pla(7),'color',COL(6,:),'linewidth',3);
            set(pla(8),'facecolor',COL(6,:));
        else
            set(pla(1:2),'color',COL(1,:),'linewidth',3);
            set(pla(3:4),'color',COL(2,:),'linewidth',3);
            set(pla(5:6),'color',COL(3,:),'linewidth',3);
            set(pla(7:8),'color',COL(6,:),'linewidth',3);
        end
    end
        
        
        if i1==i2-1&&NP<=4
            xlabel([obj.MCMC_Chain_Obj.par_names{PARINDS(i2)}]);
            ylabel([obj.MCMC_Chain_Obj.par_names{PARINDS(i1)}]);
        end
        
        for i=1:2
            YLIM = get(gca,'ylim');
            XLIM = get(gca,'xlim');
            if ~isempty(BEST)
                plot(XLIM,BEST(i1)*[1,1],'k--')
                plot(BEST(i2)*[1,1],YLIM,'k--')
            end
        end
        if NP>=4
            set(gca,'xticklabel',[],'yticklabel',[],'fontsize',9)
%             pos = get(gca,'position');
%             set(gca,'position',[pos(1:2),pos(3:4)*1.2]);
        else
            legend([pla(9),pla(7),pla(5),pla(3),pla(1)],{'means','vars','dists','vars4','Gaussian'})
        end

    end
end
if NP>=4
    for i1=1:NP-1
        for i2=i1+1:NP
            subplot(NP-1,NP-1,(i1-1)*(NP-1)+i2-1);
            pos = get(gca,'position');
            set(gca,'position',[pos(1:2)-[0.1 -0.02],pos(3:4).*[1.85 2]]);
        end
    end
end
        
        

function Make_Bias_Uncertainty_Plot(Hog,f,g)
%%
figure(f);clf
set(f,'position',[1000         627         435         711]);
set(f,'name','Parameter Uncertainty and Bias');

if ~strcmp(Hog.MOD.Gene,'STL1_and_CTT1')
    subplot(3,1,1)
    TMPmin = 10.^min(Hog.MC_Results.Bias(:,2:3),[],2);
    TMPmax = 10.^max(Hog.MC_Results.Bias(:,2:3),[],2);
    bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]); hold on;
    bar([TMPmax(1:floor(end/2)),TMPmax(ceil(end/2)+1:end)],'r');
    bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]);
    set(gca,'yscale','log','fontsize',14,'xticklabel',{'\mu','\mu,\sigma','P','4th'},'ylim',[1e-5 10^15])
    ylabel('Parameter Bias')
    legend('NonSpatial','Spatial')
   
    subplot(3,1,3)
    TMPmin = min(sqrt(Hog.MC_Results.Uncert(:,2:3).^2+(10.^Hog.MC_Results.Bias(:,2:3)).^2),[],2);
    TMPmax = max(sqrt(Hog.MC_Results.Uncert(:,2:3).^2+(10.^Hog.MC_Results.Bias(:,2:3).^2)),[],2);
    bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]); hold on;
    bar([TMPmax(1:floor(end/2)),TMPmax(ceil(end/2)+1:end)],'r');
    bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]);
    set(gca,'yscale','log','fontsize',14,'xticklabel',{'\mu','\mu,\sigma','P','4th'},'ylim',[1e-5 10^15])
    ylabel('Total Error')
    xlabel('Analysis type')
    legend('NonSpatial','Spatial')
    
end
% set(gca,'yscale','log','fontsize',14,'xticklabel',{'\mu','\mu,\sigma','P'},'ylim',[1e-5 10^15])

TMPmin = min(Hog.MC_Results.Uncert(:,2:3),[],2);
TMPmax = max(Hog.MC_Results.Uncert(:,2:3),[],2);
if ~strcmp(Hog.MOD.Gene,'STL1_and_CTT1')
    subplot(3,1,2)
    bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]); hold on;
    bar([TMPmax(1:floor(end/2)),TMPmax(ceil(end/2)+1:end)],'r');
    bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]);
    set(gca,'yscale','log','fontsize',14,'xticklabel',{'\mu','\mu,\sigma','P','4th'},'ylim',[1e-5 10^15])
else
    figure(f);clf
    set(f,'position',[228 1036 998 278]);
    set(f,'name','Parameter Uncertainty and Bias');
    bar([TMPmin([1:3,7:9,13,14]),TMPmin([4:6,10:12,15,16])]); hold on;
    bar([TMPmax([1:3,7:9,13,14]),TMPmax([4:6,10:12,15,16])],'r');
    bar([TMPmin([1:3,7:9,13,14]),TMPmin([4:6,10:12,15,16])]);
    set(gca,'yscale','log','fontsize',14,'xticklabel',{'\mu_{STL1}','\mu_{STL1},\sigma_{STL1}','P_{STL1}','\mu_{CTT1}','\mu_{CTT1},\sigma_{CTT1}','P_{CTT1}','Joint \mu','Joint \mu,\sigma'},'ylim',[1e-5 10^5])
end
ylabel('Parameter Uncertainty')
legend('NonSpatial','Spatial')


figure(f+1);clf
if ~strcmp(Hog.MOD.Gene,'STL1_and_CTT1')
    for i=1:size(Hog.MC_Results.Bias_Indv1,2)
        subplot(3,3,i)
        TMPmin = 10.^min(Hog.MC_Results.Bias_Indv1(:,i),Hog.MC_Results.Bias_Indv2(:,i));
        TMPmax = 10.^max(Hog.MC_Results.Bias_Indv1(:,i),Hog.MC_Results.Bias_Indv2(:,i));
        bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]); hold on;
        bar([TMPmax(1:floor(end/2)),TMPmax(ceil(end/2)+1:end)],'r');
        bar([TMPmin(1:floor(end/2)),TMPmin(ceil(end/2)+1:end)]);
        YLIM = get(gca,'ylim');
        YLIM(1)=1; YLIM(2) = 10^ceil(log10(YLIM(2)));
        set(gca,'yscale','log','fontsize',14,'xticklabel',[],'ylim',YLIM)
        if i==1||i==4||i==7
            ylabel('Parameter Bias')
        end
        if i==7||i==8||i==9
        set(gca,'yscale','log','fontsize',14,'xticklabel',{'\mu','\mu,\sigma','4th'},'ylim',YLIM)
        end
    end
end


