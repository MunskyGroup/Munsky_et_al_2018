classdef Hog_Model
    properties
        Parameter_File = [];
        Spatial = 'NonSpatial';
        Type = 'Moments';
        Fit_Type = 'Moments';
        Moment_Order = 2;
        Gene = 'STL1';
        MGenes = 1;
        Reps = '_AC';
        Model = '4State_k21only_noloop';
        Hog1p_Type = 'WT'; % Type of Hog1p signal {'WT','ARP8','GCN5','Hot1-5x'}
        Salt = 0.2;
        Data_File = 'Hog_Data/Analyzed_FISH_Data_2016_06';
        on_thresh = 3;
        SSA_Runs = 1000; % Number of SSA runs.
        SSA_Expts = 1; % Number of SSA experiments.
        SSA_Cells = 100; % Number of SSA cells.
        Recompile_MCMC = 0;
        MCMC_File = [];
        MCMC_Raw_Dir = [];
        gene_length = 1000; % Nt
        MCMC_Reps=1;
        k_elong = 1;  % lower limit on Nt per s
        f = 0.00;  % Correction for probability that two random spots will overlap.
        tt=[0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];  % time of data points in minutes.
        TS_Lim = 100;
        Model_Override=0;
        Mod_Ovrd_Struct=[];
    end
    properties (Constant)
        Num_States = 4;
        Search = 'no';
        Hog_Signal = 'Default';
    end
    properties (Dependent)
        Model_Obj
        Data
        MCMC_Chain_Obj
    end
    methods
        function xModel_Obj = get.Model_Obj(obj)
            if obj.Model_Override
                xModel_Obj = obj.Mod_Ovrd_Struct;
            else
                path('model_defs',path)
                xModel_Obj = Def_Model_Rules_Class(obj);
            end
        end
        function xData = get.Data(obj)
            xData = Process_Hog_Data_Class(obj);  % Call simple routine to specify the data to be fit.
        end
        function xMCMC_Chain_Obj = get.MCMC_Chain_Obj(obj)
            xMCMC_Chain_Obj = Process_MCMC_class(obj);  % Call simple routine to specify the data to be fit.
        end
    end
    methods (Static)
        function xSolutions = solve_hog(obj)
            path('solvers',path);
            path('solvers/moments',path);
            path('solvers/fsp',path);
            path('solvers/fsp',path);
            path('solvers/ssa',path);
            xSolutions = solve_hog_class(obj);
        end
        function plot_moments_mod(obj,Solutions,fignum1,fignum2)
            salts = obj.Salt;
            for isalt=1:length(obj.Salt)
                obj.Salt = salts(isalt);
                %% Means plot
                f=figure(fignum1);
                set(f,'Name','Means vs time')
                subplot(1,length(salts),isalt)
                switch obj.Spatial
                    case {'NonSpatial'}
                        plot(obj.tt,max(1e-6,Solutions(isalt).moments.Trajectories(:,obj.Num_States+1:obj.Num_States+obj.MGenes))); hold on
                    case {'Spatial','Fixed'}
                        plot(obj.tt,max(1e-6,Solutions(isalt).moments.Trajectories(:,obj.Num_States+1:obj.Num_States+2*obj.MGenes))); hold on
                end
                set(get(gca,'children'),'linewidth',3);
                set(gca,'fontsize',14,'yscale','log')
                xlabel('time (min)')
                ylabel('means');
                %% Vars plots
                if obj.Moment_Order>=2
                    f=figure(fignum2);
                    set(f,'Name','VARs vs time')
                    subplot(1,length(salts),isalt)
                    switch obj.Spatial
                        case {'NonSpatial'}
                            inds_vars = obj.Num_States+obj.MGenes+... %mns
                                sum((obj.Num_States+obj.MGenes):-1:obj.MGenes+1)+... %gene states
                                (1:sum(1:obj.MGenes));
                            %                             plot(obj.tt,Solutions(isalt).moments.Trajectories(:,inds_vars));hold on;
                            plot(obj.tt,max(1e-6,sqrt(Solutions(isalt).moments.Trajectories(:,inds_vars))));hold on;
                        case {'Spatial','Fixed'}
                            inds_vars = obj.Num_States+obj.MGenes*2+... %mns
                                sum(obj.Num_States+2*obj.MGenes:-1:2*obj.MGenes+1)+... %gene states
                                (1:sum(1:obj.MGenes*2));
                            %                             plot(obj.tt,Solutions(isalt).moments.Trajectories(:,inds_vars));hold on;
                            plot(obj.tt,max(1e-6,sqrt(Solutions(isalt).moments.Trajectories(:,inds_vars))));hold on;
                    end
                    set(get(gca,'children'),'linewidth',3);
                    set(gca,'fontsize',14,'yscale','log')
                    xlabel('time (min)')
                    ylabel('(co)variances');
                end
            end
        end
        function plot_moments_dat(obj,Data,fignum1,fignum2)
            salts = obj.Salt;
            for isalt=1:length(obj.Salt)
                %                 obj.Salt = salts(isalt);
                %% Means plot
                figure(fignum1);
                subplot(1,length(salts),isalt)
                switch obj.Spatial
                    case {'NonSpatial'}
                        plot(Data(isalt).Times/60,Data(isalt).Trajectories(:,1:obj.MGenes),'-o');hold on
                    case {'Spatial','Fixed'}
                        plot(Data(isalt).Times/60,Data(isalt).Trajectories(:,1:2*obj.MGenes),'-o');hold on
                end
                %                 set(get(gca,'children'),'linewidth',3);
                set(gca,'fontsize',14,'yscale','log')
                title(['Means, Salt = ',num2str(salts(isalt))])
                xlabel('time (min)')
                ylabel('means');
                %% Vars plots
                if obj.Moment_Order>=2
                    figure(fignum2);
                    subplot(1,length(salts),isalt)
                    switch obj.Spatial
                        case {'NonSpatial'}
                            %                             plot(Data(isalt).Times/60,Data(isalt).Trajectories(:,obj.MGenes+1:1:end),'-o');hold on
                            plot(Data(isalt).Times/60,sqrt(Data(isalt).Trajectories(:,obj.MGenes+1:1:end)),'-o');hold on
                        case {'Spatial','Fixed'}
                            %                             plot(Data(isalt).Times/60,Data(isalt).Trajectories(:,2*obj.MGenes+1:end),'-o');hold on
                            plot(Data(isalt).Times/60,sqrt(Data(isalt).Trajectories(:,2*obj.MGenes+1:end)),'-o');hold on
                    end
                    title(['(co)Variances, Salt = ',num2str(salts(isalt))])
                    xlabel('time (min)')
                    ylabel('(co)variances');
                    %                     set(get(gca,'children'),'linewidth',3);
                    set(gca,'fontsize',14,'yscale','log')
                end
            end
        end
        function plot_distributions_mod(obj,Solutions,fignum1,fignum2)
            for isalt=1:length(obj.Salt)
                f = figure(fignum1-1+isalt);
                set(f,'Name',[obj.Gene,' Dists, ',num2str(obj.Salt(isalt))]);
                switch(obj.Spatial)
                    case 'NonSpatial'
                        Y_MOD = squeeze(sum(Solutions(isalt).distributions,3));
                    case {'Spatial','Fixed'}
                        Y_MOD = zeros(size(Solutions(isalt).distributions,1),...
                            size(Solutions(isalt).distributions,2)+...
                            size(Solutions(isalt).distributions,3)-1);
                        for i=1:size(Solutions(isalt).distributions,1)
                            for j=1:size(Solutions(isalt).distributions,2)
                                for k=1:size(Solutions(isalt).distributions,3)
                                    Y_MOD(i,j+k-1) = Y_MOD(i,j+k-1)+Solutions(isalt).distributions(i,j,k);
                                end
                            end
                        end
                end
                on_v_t = zeros(size(obj.tt));
                for i_time = 1:length(obj.tt)
                    subplot(4,4,i_time);
                    TMP = size(Y_MOD,2);
                    plot(0:TMP-1,smooth_hists(Y_MOD(i_time,:),3),'linewidth',3); hold on;
                    set(gca,'xlim',[0 150],'ylim',[0 0.04],'yticklabel',[],'xticklabel',[]);
                    on_v_t(i_time) = sum(Y_MOD(i_time,obj.on_thresh+1:end));
                end
                f2=figure(fignum2);
                set(f2,'Name',[obj.Gene,' on vs. time, ',num2str(obj.Salt(isalt))]);
                subplot(1,2,isalt);
                plot(obj.tt,max(1e-6,on_v_t),'linewidth',3); hold on
                set(gca,'xlim',[0 60],'ylim',[0 1],'fontsize',16);
            end
           
            for isalt=1:length(obj.Salt)
%                 fa = figure(1234+isalt);
                switch(obj.Spatial)
                    case 'NonSpatial'
                        Y_MOD = squeeze(sum(Solutions(isalt).distributions,3));
                    case {'Spatial','Fixed'}
                        Y_MOD = zeros(size(Solutions(isalt).distributions,1),...
                            size(Solutions(isalt).distributions,2)+...
                            size(Solutions(isalt).distributions,3)-1);
                        for i=1:size(Solutions(isalt).distributions,1)
                            for j=1:size(Solutions(isalt).distributions,2)
                                for k=1:size(Solutions(isalt).distributions,3)
                                    Y_MOD(i,j+k-1) = Y_MOD(i,j+k-1)+Solutions(isalt).distributions(i,j,k);
                                end
                            end
                        end
                end
%                 figure(fa)
%                 Q = [1,7,9,11,13];
%                 for v = 1:length(Q)%1:length(obj.tt)
%                     i_time = Q(v);
%                     TMP = size(Y_MOD,2);
%                     plot3(obj.tt(i_time)*ones(1,TMP),0:TMP-1,smooth_hists(Y_MOD(i_time,:),3),'k','linewidth',3); hold on;
%                     set(gca,'xlim',[0 40],'ylim',[0 100],'zlim',[0 0.04],'yticklabel',[],'xticklabel',[]);
%                     set(gca,'xtick',obj.tt(Q),'xticklabel',obj.tt(Q),'xdir','reverse')
%                 end
                
                fb = figure(1236+isalt);
                figure(fb)
                set(fb,'position',[514   754   578    91]);
                Q = [1,7,9,11,13];
                for v = 1:length(Q)%1:length(obj.tt)
                    subplot(1,5,v)
                    i_time = Q(v);
                    TMP = size(Y_MOD,2);
                    plot(0:TMP-1,smooth_hists(Y_MOD(i_time,:),3),'k','linewidth',3); hold on;
%                     set(gca','position',[0.05+0.13*(v-1) 0.1100 0.1237 0.8150],'yticklabel',[],'xticklabel',[])
                    set(gca,'xlim',[0 100],'ylim',[0 0.04],'yticklabel',[],'xticklabel',[]);
%                     set(gca,'xtick',obj.tt(Q),'xticklabel',obj.tt(Q),'xdir','reverse')
                end
            end

            
        end
        function plot_distributions_dat(obj,Data,fignum1,fignum2)
            for isalt=1:length(obj.Salt)
                switch(obj.Spatial)
                    case 'NonSpatial'
                        Y_DAT = Data(isalt).Trajectories;
                    case {'Spatial','Fixed'}
                        Y_DAT = zeros(size(Data(isalt).Trajectories,1),...
                            size(Data(isalt).Trajectories,2)+...
                            size(Data(isalt).Trajectories,3)-1);
                        for i=1:size(Data(isalt).Trajectories,1)
                            for j=1:size(Data(isalt).Trajectories,2)
                                for k=1:size(Data(isalt).Trajectories,3)
                                    Y_DAT(i,j+k-1) = Y_DAT(i,j+k-1)+Data(isalt).Trajectories(i,j,k);
                                end
                            end
                        end
                end
                on_v_t = zeros(size(obj.tt));
                figure(fignum1-1+isalt);
                for i_time = 1:length(obj.tt)
                    subplot(4,4,i_time);
                    TMP = size(Y_DAT,2);
                    nrm = max(1,Data(isalt).Num_cells(i_time));
                    plot(0:TMP-1,smooth_hists(Y_DAT(i_time,:)/nrm,3),'linewidth',1); hold on;
                    set(gca,'xlim',[0 150],'ylim',[0 0.04],'yticklabel',[],'xticklabel',[]);
                    on_v_t(i_time) = sum(Y_DAT(i_time,obj.on_thresh+1:end)/nrm);
                end
                figure(fignum2);
                subplot(1,2,isalt);
                plot(obj.tt,on_v_t,'o');
            end
            
            for isalt=1:length(obj.Salt)
                switch(obj.Spatial)
                    case 'NonSpatial'
                        Y_DAT = Data(isalt).Trajectories;
                    case {'Spatial','Fixed'}
                        Y_DAT = zeros(size(Data(isalt).Trajectories,1),...
                            size(Data(isalt).Trajectories,2)+...
                            size(Data(isalt).Trajectories,3)-1);
                        for i=1:size(Data(isalt).Trajectories,1)
                            for j=1:size(Data(isalt).Trajectories,2)
                                for k=1:size(Data(isalt).Trajectories,3)
                                    Y_DAT(i,j+k-1) = Y_DAT(i,j+k-1)+Data(isalt).Trajectories(i,j,k);
                                end
                            end
                        end
                end
%                 fa = figure(1234+isalt);
%                 figure(fa)
%                 set(gcf,'Position',[1956         766         211         170]);
%                 Q = [1,3,4,5,6,7,8,9,10,11,12,13];
                 Q = [1,7,9,11,13];
%                 for v = 1:length(Q)%1:length(obj.tt)
%                     i_time = Q(v);
%                     TMP = size(Y_DAT,2);
%                     nrm = max(1,Data(isalt).Num_cells(i_time));
%                     plot3(obj.tt(i_time)*ones(1,TMP),0:TMP-1,smooth_hists(Y_DAT(i_time,:)/nrm,3),'r','linewidth',2); hold on;
%                     set(gca,'xlim',[0 40],'ylim',[0 100],'zlim',[0 0.04],'yticklabel',[],'xticklabel',[]);
%                 end
%                 set(gca,'Position',[ 0.1885    0.2131    0.6613    0.7119],...
%                     'view',[ 164.1000   23.6000]);
%                 set(gca,'xtick',obj.tt(Q),'xticklabel',obj.tt(Q),'xdir','reverse')
%                 set(gca,'ytick',[0 50 100 150])
%                 set(gca,'ytick',[0 50 100 150],'yticklabel',[0 50 100 150])
%                 set(gca,'box','on')

                fb = figure(1236+isalt); figure(fb)
                set(fb,'position',[514   754   578    91]);
                for v = 1:length(Q)%1:length(obj.tt)
                    subplot(1,5,v)
                    i_time = Q(v);
                    TMP = size(Y_DAT,2);
                    nrm = max(1,Data(isalt).Num_cells(i_time));
                    plot(0:TMP-1,smooth_hists(Y_DAT(i_time,:)/nrm,3),'r','linewidth',2); hold on;
                    set(gca,'xlim',[0 100],'ylim',[0 0.04],'yticklabel',[],'xticklabel',[]);
                    %                     set(gca,'xtick',obj.tt(Q),'xticklabel',obj.tt(Q),'xdir','reverse')
                end
            end
        end
        function plot_spat_distributions_mod(obj,Solutions,fignum1,fignum2)
            bin = 10;
            Thresh = 0.01;
            Thresh_Low = 0.001;
            for isalt=1:length(obj.Salt)
                f = figure(fignum1-1+isalt);
                set(f,'Name',[obj.Gene,' Spat. Dists, ',num2str(obj.Salt(isalt))]);
                
                Y_MOD = Solutions(isalt).distributions;
                
                Nb3 = min(size(Y_MOD,3),150);
                Nb2 = min(size(Y_MOD,2),150);
                %                 Chunk2 = min(bin,max(1,floor(Nb2/15)))
                %                 Chunk3 = min(bin,max(1,floor(Nb3/15)))
                Chunk2 = 1;
                Chunk3 = bin;
                on_v_t = zeros(size(Y_MOD,1),3);
                for i_time = 1:size(Y_MOD,1)
                    subplot(4,4,i_time);
                    on_v_t(i_time,1) = 1-sum(sum(squeeze(Y_MOD(i_time,1:obj.on_thresh,:))));
                    on_v_t(i_time,2) = 1-sum(sum(squeeze(Y_MOD(i_time,:,1:obj.on_thresh))));
                    on_v_t(i_time,3) = sum(sum(squeeze(Y_MOD(i_time,obj.on_thresh+1:end,obj.on_thresh+1:end))));
                    
                    TMP = squeeze(Y_MOD(i_time,:,:));
                    i=1;
                    clear TMP2;
                    while i*Chunk3<size(TMP,2);
                        j=1;
                        while j*Chunk2<size(TMP,1);
                            TMP2(j,i) = sum(sum(TMP((j-1)*Chunk2+1:j*Chunk2,(i-1)*Chunk3+1:i*Chunk3),2));
                            j=j+1;
                        end
                        i=i+1;
                    end
                    contourf((0:size(TMP2,2)-1)*Chunk3,(0:size(TMP2,1)-1)*Chunk2,-max(Thresh_Low,min(Thresh,TMP2)))
                    colormap('gray')
                end
                f2=figure(fignum2);
                subplot(1,2,isalt);
                set(f2,'Name',['spatial on vs. time, ',num2str(obj.Salt(isalt))]);
                plot(obj.tt,on_v_t,'linewidth',3); hold on
                set(gca,'xlim',[0 60],'ylim',[0 1],'fontsize',16);
            end
        end
        function plot_spat_distributions_dat(obj,Data,fignum1,fignum2)
            bin = 10;
            Thresh = 0.01;
            Thresh_Low = 0.001;
            for isalt=1:length(obj.Salt)
                f = figure(fignum1-1+isalt);
                
                Y_DAT = Data(isalt).Trajectories;
                %                 Nb3 = min(size(Y_DAT,3),150);
                %                 Nb2 = min(size(Y_DAT,2),150);
                %                 Chunk2 = min(bin,max(1,floor(Nb2/15)))
                %                 Chunk3 = min(bin,max(1,floor(Nb3/15)))
                %                 pause
                Chunk2 = 1;
                Chunk3 = bin;
                on_v_t = zeros(size(Y_DAT,1),3);
                for i_time = 1:size(Y_DAT,1)
                    nrm = max(1,Data(isalt).Num_cells(i_time));
                    subplot(4,4,i_time);
                    if sum(sum(squeeze(Y_DAT(i_time,:,:))))>0
                        on_v_t(i_time,1) = 1-sum(sum(squeeze(Y_DAT(i_time,1:obj.on_thresh,:))))/nrm;
                        on_v_t(i_time,2) = 1-sum(sum(squeeze(Y_DAT(i_time,:,1:obj.on_thresh))))/nrm;
                        on_v_t(i_time,3) = sum(sum(squeeze(Y_DAT(i_time,obj.on_thresh+1:end,obj.on_thresh+1:end))))/nrm;
                    end
                    TMP = squeeze(Y_DAT(i_time,:,:))/nrm;
                    i=1;
                    TMP2 = Thresh_Low*ones(15,20);
                    while i*Chunk3<size(TMP,2);
                        j=1;
                        while j*Chunk2<size(TMP,1);
                            TMP2(j,i) = sum(sum(TMP((j-1)*Chunk2+1:j*Chunk2,(i-1)*Chunk3+1:i*Chunk3),2));
                            j=j+1;
                        end
                        i=i+1;
                    end
                    contourf((0:size(TMP2,2)-1)*Chunk3,(0:size(TMP2,1)-1)*Chunk2,-max(Thresh_Low,min(Thresh,TMP2)))
                    colormap('gray')
                end
                figure(fignum2);
                subplot(1,2,isalt);
                plot(obj.tt,on_v_t,'o'); hold on
                set(gca,'xlim',[0 60],'ylim',[0 1],'fontsize',16);
            end
        end
        function plot_SSA_moment_results(obj,SSA,fignum1,fignum2)
            Nt = length(obj.tt);
            %             TMPm = zeros(obj.SSA_Expts,length(obj.Salt),Nt,);
            for ks=1:obj.SSA_Expts
                for isalt=1:length(obj.Salt)
                    SSA_RNA = zeros(obj.SSA_Cells,size(SSA(isalt).SSA,2),Nt);
                    for it=1:Nt
                        SMPs = ceil(obj.SSA_Runs*rand(obj.SSA_Cells,1));
                        SSA_RNA(:,:,it)=squeeze(SSA(isalt).SSA(SMPs,:,it));
                        %% Corrrection for spot overlaps.
                        for j=1:size(SSA_RNA,1)
                            for j2=1:size(SSA_RNA,2)
                                c=0;
                                kl = SSA_RNA(j,j2,it);
                                for z=1:kl
                                    if rand<(1-obj.f)^c
                                        c=c+1;
                                    end
                                end
                                SSA_RNA(j,j2,it)=c;
                            end
                        end
                        MN = mean(squeeze(SSA_RNA(:,:,it)),1);
                        COV = cov(squeeze(SSA_RNA(:,:,it)));
                        Solutions(isalt).moments.Trajectories(it,1:length(MN)) = MN;
                        k = length(MN)+1;
                        for i=1:length(MN)
                            for j=i:length(MN)
                                Solutions(isalt).moments.Trajectories(it,k) = COV(i,j);
                                k=k+1;
                            end
                        end
                    end
                    TMPm(ks,isalt,:,:) = Solutions(isalt).moments.Trajectories;
                    obj.plot_moments_mod(obj,Solutions,fignum1,fignum2);
                end
            end
            Solutions_Med=Solutions;
            Solutions_Med.moments.Trajectories = squeeze(median(TMPm(:,1,:,:)));
            obj.plot_moments_mod(obj,Solutions_Med,fignum1,fignum2);
        end
        function plot_SSA_dist_results(obj,SSA,fignum1,fignum2,sp1,sp2)
            Nt = length(obj.tt);
            for ks=1:obj.SSA_Expts
                for isalt=1:length(obj.Salt)
                    SSA_RNA = zeros(obj.SSA_Cells,size(SSA(isalt).SSA,2),Nt);
                    for it=1:Nt
                        SMPs = ceil(obj.SSA_Runs*rand(obj.SSA_Cells,1));
                        SSA_RNA(:,:,it)=squeeze(SSA(isalt).SSA(SMPs,:,it));
                        %% Corrrection for spot overlaps.
                        for j=1:size(SSA_RNA,1)
                            for j2=1:size(SSA_RNA,2)
                                c=0;
                                kl = SSA_RNA(j,j2,it);
                                for z=1:kl
                                    if rand<(1-obj.f)^c
                                        c=c+1;
                                    end
                                end
                                SSA_RNA(j,j2,it)=c;
                            end
                        end
                    end
                    Nrm = size(SSA_RNA,1);
                    if exist('sp2','var')&&~isempty(sp2)
                        mx1 = max(max(squeeze(sum(SSA_RNA(:,obj.Num_States+sp1,:),2))));
                        mx2 = max(max(squeeze(sum(SSA_RNA(:,obj.Num_States+sp2,:),2))));
                        Solutions(isalt).distributions = zeros(Nt,mx1+1,mx2+1);
                        for it = 1:Nt
                            for i=1:mx1
                                for j=1:mx2
                                    Solutions(isalt).distributions(it,i,j) = ...
                                        sum((sum(SSA_RNA(:,obj.Num_States+sp1,it),2)==i-1)&...
                                        (sum(SSA_RNA(:,obj.Num_States+sp2,it),2)==j-1))/Nrm;
                                end
                            end
                        end
                    else
                        mx1 = max(max(squeeze(sum(SSA_RNA(:,obj.Num_States+sp1,:),2))));
                        Solutions(isalt).distributions = zeros(Nt,mx1+1,1);
                        for it = 1:Nt
                            for i=1:mx1
                                Solutions(isalt).distributions(it,i,1) = ...
                                    sum(sum(SSA_RNA(:,obj.Num_States+sp1,it),2)==i-1)/Nrm;
                            end
                        end
                    end
                    TMPm(isalt).distributions(ks,1:Nt,1:mx1+1,1) = Solutions(isalt).distributions;
                end
                
                if exist('sp2','var')&&~isempty(sp2)
                    obj.plot_spat_distributions_mod(obj,Solutions,fignum1,fignum2)
                else
                    obj.plot_distributions_mod(obj,Solutions,fignum1,fignum2)
                end
            end
            Solutions_Med=Solutions;
            for ksalt = 1:isalt
                Solutions_Med(ksalt).distributions = squeeze(median(TMPm(isalt).distributions));
                if exist('sp2','var')&&~isempty(sp2)
                    obj.plot_spat_distributions_mod(obj,Solutions_Med,fignum1,fignum2)
                else
                    obj.plot_distributions_mod(obj,Solutions_Med,fignum1,fignum2)
                end
            end
        end
        function MCMC_Results = plot_MCMC_results(obj,MCMC,show_pars,fignum1)
            path('solvers',path);
            
            
            if length(MCMC.chain)==2
                A_Chain.pars = MCMC.chain(1).pars;
                B_Chain.pars = MCMC.chain(2).pars;
            else
                A_Chain.pars=[];
                for i=1:floor(length(MCMC.chain)/2)
                    A_Chain.pars = [A_Chain.pars;MCMC.chain(i).pars];
                end
                B_Chain.pars=[];
                for i=floor(length(MCMC.chain)/2)+1:length(MCMC.chain)
                    B_Chain.pars = [B_Chain.pars;MCMC.chain(i).pars];
                end
            end
            
            COV = cov(log10(abs([A_Chain.pars;...
                B_Chain.pars])),1);
            LogM = mean(log10(abs([A_Chain.pars;...
                B_Chain.pars])),1);
            
            NP = length(show_pars);
            f = figure(fignum1);
            set(f,'Name',[obj.Gene,', ',obj.Spatial,', MCMC Scatter']);
            COV1 = cov(log10([A_Chain.pars]));
            COV2 = cov(log10([B_Chain.pars]));
            LogM1 = mean(log10(abs([A_Chain.pars])),1);
            LogM2 = mean(log10(abs([B_Chain.pars])),1);
            
            if obj.MCMC_Reps
                for i1=1:NP-1
                    for i2=i1+1:NP
                        subplot(NP-1,NP-1,(i1-1)*(NP-1)+i2-1);
                        cv = [COV1(show_pars(i2),show_pars(i2)),COV1(show_pars(i2),show_pars(i1));...
                            COV1(show_pars(i1),show_pars(i2)),COV1(show_pars(i1),show_pars(i1))];
                        mu = [LogM1(show_pars(i2));LogM1(show_pars(i1))];
                        error_ellipse_empty(cv,mu,'conf',0.9)  %% Call routine to plots ellipses for the 90% confidense intervals.
                        hold on
                        
                        subplot(NP-1,NP-1,(i1-1)*(NP-1)+i2-1);
                        cv = [COV2(show_pars(i2),show_pars(i2)),COV2(show_pars(i2),show_pars(i1));...
                            COV2(show_pars(i1),show_pars(i2)),COV2(show_pars(i1),show_pars(i1))];
                        mu = [LogM2(show_pars(i2));LogM2(show_pars(i1))];
                        error_ellipse_empty(cv,mu,'conf',0.9)  %% Call routine to plots ellipses for the 90% confidense intervals.
                    end
                end
            else
                for i1=1:NP-1
                    for i2=i1+1:NP
                        subplot(NP-1,NP-1,(i1-1)*(NP-1)+i2-1);
                        cv = [COV(show_pars(i2),show_pars(i2)),COV(show_pars(i2),show_pars(i1));...
                            COV(show_pars(i1),show_pars(i2)),COV(show_pars(i1),show_pars(i1))];
                        mu = [LogM(show_pars(i2));LogM(show_pars(i1))];
                        error_ellipse_fill(cv,mu,'conf',0.9)  %% Call routine to plots ellipses for the 90% confidense intervals.
                        hold on
                    end
                end
            end
            
            f2 = figure(fignum1+1);
            set(f2,'Name',[obj.Gene,', ',obj.Spatial,', MCMC Scatter']);
            COV1 = cov(log10([A_Chain.pars]));
            COV2 = cov(log10([B_Chain.pars]));
            for i1=1:NP-1
                for i2=i1+1:NP
                    subplot(NP-1,NP-1,(i1-1)*(NP-1)+i2-1);
                    cv = [COV1(show_pars(i2),show_pars(i2)),COV1(show_pars(i2),show_pars(i1));...
                        COV1(show_pars(i1),show_pars(i2)),COV1(show_pars(i1),show_pars(i1))];
                    mu = [0;0];
                    error_ellipse_empty(cv,mu,'conf',0.9)  %% Call routine to plots ellipses for the 90% confidense intervals.
                    hold on
                    
                    cv = [COV2(show_pars(i2),show_pars(i2)),COV2(show_pars(i2),show_pars(i1));...
                        COV2(show_pars(i1),show_pars(i2)),COV2(show_pars(i1),show_pars(i1))];
                    mu = [0;0];
                    error_ellipse_empty(cv,mu,'conf',0.9)  %% Call routine to plots ellipses for the 90% confidense intervals.
                    hold on
                end
            end

            
            if nargout==1
                COV = cov([A_Chain.pars;...
                    B_Chain.pars]);
                COV1 = cov([A_Chain.pars]);
                COV2 = cov([B_Chain.pars]);
                
                %                 COV = cov([A_Chain.pars;...
                %                     B_Chain.pars]./...
                %                     repmat(abs(mean([A_Chain.pars;...
                %                     B_Chain.pars],1)),...
                %                     size(A_Chain.pars,1)+size(B_Chain.pars,1),1),1);
                %                 COV1 = cov([A_Chain.pars]./...
                %                     repmat(abs(mean([A_Chain.pars;...
                %                     B_Chain.pars],1)),size(A_Chain.pars,1),1),1);
                %                 COV2 = cov([B_Chain.pars]./...
                %                     repmat(abs(mean([A_Chain.pars;...
                %                     B_Chain.pars],1)),size(B_Chain.pars,1),1),1);

                
%                 COV_log = cov(log10(abs([A_Chain.pars;...
%                 B_Chain.pars])),1);
                COV_log = log10(abs(cov([A_Chain.pars;...
                B_Chain.pars],1)));
                COV_sign = sign(COV);

                MCMC_Results.COV_log = COV_log(show_pars,show_pars);
                MCMC_Results.COV_sign = COV_sign(show_pars,show_pars);
                
                MCMC_Results.trace = trace(COV(show_pars,show_pars));
                MCMC_Results.trace1 = trace(COV1(show_pars,show_pars));
                MCMC_Results.trace2 = trace(COV2(show_pars,show_pars));
                
                MCMC_Results.M1 = mean((abs(A_Chain.pars(:,show_pars))),1);
                MCMC_Results.M2 = mean((abs(B_Chain.pars(:,show_pars))),1);
                MCMC_Results.M = mean((abs([A_Chain.pars(:,show_pars);...
                    B_Chain.pars(show_pars)])),1);
                
                DD = [A_Chain.pars(:,show_pars);B_Chain.pars(show_pars)];
                MCMC_Results.STD = std((abs(DD)),0,1);
                
                L = size(DD,1);
                for ij = 1:length(show_pars)
                    a = sort(DD(:,ij));
                    MCMC_Results.LB99(ij) = a(floor(0.01*L));
                    MCMC_Results.UB99(ij) = a(ceil(0.99*L));
                    MCMC_Results.LB90(ij) = a(floor(0.05*L));
                    MCMC_Results.UB90(ij) = a(ceil(0.95*L));
                    MCMC_Results.LB50(ij) = a(floor(0.25*L));
                    MCMC_Results.UB50(ij) = a(ceil(0.75*L));
                end
            end
            
        end
    end
end
