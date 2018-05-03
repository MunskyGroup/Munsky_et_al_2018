clear all
close all
%% MCMC Figures
MOD = Hog_Model;
MOD.Model = '4State_k21only_noloop';
MOD.Recompile_MCMC = 1;
%%
col = ['r';'b';'k';'m'];
for i = 1:2
    switch i
        case 1
            MOD.Gene = 'STL1';
        case 2
            MOD.Gene = 'CTT1';
    end
    for k=1:2
        switch k
            case 1
                MOD.Spatial = 'NonSpatial';
            case 2
                MOD.Spatial = 'Spatial';
        end

        for j=1:4
            switch j
                case 1
                    MOD.Fit_Type = 'Means';
                case 2
                    MOD.Fit_Type = 'Moments_Corrected';
                case 3
                    MOD.Fit_Type = 'Distributions';
                case 4
                    MOD.Fit_Type = 'Moments_4';
            end

            MOD.MCMC_Raw_Dir = ['/Users/munsky/Dropbox/Shared_Codes/SSIT_Data_and_Results/Hog_System/MCMC_Results/',MOD.Gene,'_0p2and0p4/',MOD.Spatial,'/',MOD.Fit_Type,'/',MOD.Model];
            DIR = ['Fit_Results/',MOD.Gene,'_0p2and0p4/',MOD.Spatial,'/',MOD.Fit_Type,'/',MOD.Model];
            MOD.MCMC_File = [DIR,'/MCMC_Processed.mat'];

            figure(10*(i-1)+1);
            subplot(4,2,(j-1)*2+k)
            TMP = MOD.MCMC_Chain_Obj;
            kj = 1;
            
            A = TMP.chain(1).obj;
            a = TMP.chain(1).pars;
            
            B = TMP.chain(2).obj;
            b = TMP.chain(2).pars;
            kj=2;
            while isempty(b)
                kj=kj+1;
                B = TMP.chain(kj).obj;
                b = TMP.chain(kj).pars;
            end
            
            
            mx = max([A;B]);
            [~,N] = hist(A-mx,(-20:0));
            Z= zeros(length(N),length(TMP.chain));
            for ii=1:length(TMP.chain)
                B = TMP.chain(ii).obj;
                [Y] = hist(B-mx,N);
                Z(:,ii) = Y'/sum(Y);
                bar(N,[Z])
            end
            
%             KS = max(abs(cumsum(Z(:,1))-cumsum(Z(:,2))));
%             text(-14,0.22,['KS = ',num2str(KS)])
            
            set(gca,'fontsize',14,'xlim',[-15 0.5],'ylim',[0 0.25]);
            xlabel('log_{10} (L/L^*)')
%             
%             f100=figure(100); 
%             set(f100,'position',[802   430   459   387]);
%             subplot(2,2,(k-1)*2+i)
%             mx = max([A;B]);
%             [X,N] = hist(A-mx,(-20:0));
%             Z= zeros(length(N),2);
%             for ii=1:2
%                 B = TMP.chain(ii).obj;
%                 [Y] = hist(B-mx,N);
%                  [i,j,k,sum(Y)*100]
%                Z(:,ii) = Y'/sum(Y);
%                 plot(N,Z(:,ii),col(j),'linewidth',3); hold on
%             end
            
%%
            f100=figure(100);
            set(f100,'position',[802   430   459   387]);
            subplot(2,2,(k-1)*2+i)
            mx = max([A;B]);
            [X,N] = hist(A-mx,(-20:0));
            Z= zeros(length(N),2);
            AA1=[];AA2=[];
            for ii=1:length(TMP.chain)
                if mod(ii,2)==0
                    AA1=[AA1;TMP.chain(ii).obj];
                else
                    AA2=[AA2;TMP.chain(ii).obj];
                end
            end
            Y = hist(AA1-mx,N);
            [i,j,k,sum(Y)*100]
            bbb = Y'/sum(Y);
            plot(N,bbb,col(j),'linewidth',3); hold on
            Y = hist(AA2-mx,N);
            [i,j,k,sum(Y)*100]
            bbb = Y'/sum(Y);
            plot(N,bbb,col(j),'linewidth',3); hold on

            %%
            
            set(gca,'fontsize',14,'xlim',[-15 0.5],'ylim',[0 0.21],'ytick',[0 0.05 0.1 0.15 0.2]);
            if k==2
                xlabel('log_{10} (L/L^*)')
            end
           
            if j==3&&i==2
                f101=figure(101);
                subplot(1,2,k);
                if k==1
                    plot(a(:,10),a(:,12),'ro'); hold on
                    plot(b(:,10),b(:,12),'bo'); hold on
                else
                    plot(a(:,10),a(:,14),'ko'); hold on
                    plot(b(:,10),b(:,14),'co'); hold on
                end
%                 set(gca,'xscale','log','yscale','log');
            end
            
            figure(10*(i-1)+2);
            subplot(4,2,(j-1)*2+k)
            for ii=1:length(TMP.chain)
                B = TMP.chain(ii).obj;
                plot(-B); hold on
            end

            figure(10*(i-1)+2+j+20*(k-1));
            for ii=1:length(TMP.par_names)
                try
                    subplot(4,4,ii)
                    [X,N] = hist(a(:,ii),20);
                    [Y] = hist(b(:,ii),N);
                    title(TMP.par_names(ii))
                    plot(N,[X',Y'],'linewidth',3)
                catch
                end
            end
            
            drawnow
        end
    end
end
