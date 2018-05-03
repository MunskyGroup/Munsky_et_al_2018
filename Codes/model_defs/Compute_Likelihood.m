function [Objective]  = Compute_Likelihood(obj,Dynamics,SIG,N_Space)
if nargin<3
    SIG=0;
end
Data = obj.Data;
PARS = obj.Model_Obj.PARS;
for i_Salt = 1:length(Data)
    PARS.Salt = obj.Salt(i_Salt);
    if strcmp(obj.Type,'Moments')
%         Dynamics(i_Salt) = Run_Hog_Model_v2(PARS);
        switch  obj.Moment_Order% Order of the moments (=1 or 2).
            case 1  %% Using the model means and the data variance
                Mu_M = Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Mu.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                Mu_D = Data(i_Salt).Trajectories.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                %                 Objective(i_Salt) = sum(sum((Mu_D-Mu_M).^2./(2*Data(i_Salt).Variance.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2)))));
                %                 MNS = (Mu_D(2:end,:)-Mu_M(2:end,:)).^2./(2*Data(i_Salt).Variance(2:end,:).*repmat(Data(i_Salt).N_Cells(2:end),1,size(Data(i_Salt).Output_Mu,2)))
                VART = Data(i_Salt).Variance(2:end,:);
                Xi2 = (Mu_D(2:end,:)-Mu_M(2:end,:)).^2./(2*VART.*repmat(Data(i_Salt).N_Cells(2:end),1,size(Data(i_Salt).Output_Mu,2)));
                Xi2 = Xi2(Data(i_Salt).N_Cells(2:end)~=0,:);
                Objective(i_Salt) = sum(sum(Xi2));
            case 2
                if SIG==1  %% Using the model means and the model variance, but only computing the mean likelihoods.
                    Mu_M = Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Mu.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                    Mu_D = Data(i_Salt).Trajectories(:,1:size(Mu_M,2)).*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                    VART = max(1e-6,Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Sig.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Sig,2)));
                    VART = VART(2:end,:);
                    if size(VART,2)==1
                        Xi2 = (Mu_D(2:end,:)-Mu_M(2:end,:)).^2./(2*VART);
                    elseif size(VART,2)==3
                        Xi2 = zeros(size(Mu_D(2:end,:)));
                        for i=1:size(Xi2,1)
                            Xi2(i) = (Mu_D(i+1,:)-Mu_M(i+1,:))*...
                                (2*[VART(i,1),VART(i,2);VART(i,2),VART(i,3)])^-1*...
                                (Mu_D(i+1,:)-Mu_M(i+1,:))'...
                                + 1/2 * log(det([VART(i,1),VART(i,2);VART(i,2),VART(i,3)]));
                        end
                    end
                    Objective_Mu = sum(sum(Xi2));
                    Objective(i_Salt) = Objective_Mu;
                elseif SIG==2  %% Using the model means and the sample variance.
                    Mu_M = Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Mu.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                    Mu_D = Data(i_Salt).Trajectories(:,1:size(Mu_M,2)).*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                    VART = Data(i_Salt).Variance(2:end,:);
                    Xi2 = (Mu_D(2:end,:)-Mu_M(2:end,:)).^2./(2*VART.*repmat(Data(i_Salt).N_Cells(2:end),1,size(Data(i_Salt).Output_Mu,2)));
                    Xi2 = Xi2(Data(i_Salt).N_Cells(2:end)~=0,:);
                    Objective(i_Salt) = sum(sum(Xi2));
                elseif SIG==0
                    Mu_M = Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Mu.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                    Mu_D = Data(i_Salt).Trajectories(:,1:size(Mu_M,2)).*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                    VART = max(1e-6,Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Sig.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Sig,2)));
                    VART = VART(2:end,:);
                    if size(VART,2)==1
%                         Xi2 = (Mu_D(2:end,:)-Mu_M(2:end,:)).^2./(2*VART);
                        Xi2 = (Mu_D(2:end,:)-Mu_M(2:end,:)).^2./(2*VART) + 1/2*log(VART);
                    elseif size(VART,2)==3
                        Xi2 = zeros(size(Mu_D(2:end,:)));
                        for i=1:size(Xi2,1)
                            %                         Xi2(i) = (Mu_D(i+1,:)-Mu_M(i+1,:))*...
                            %                             (2*[VART(i,1),VART(i,2);VART(i,2),VART(i,3)])^-1*...
                            %                             (Mu_D(i+1,:)-Mu_M(i+1,:))';
                            Xi2(i) = (Mu_D(i+1,:)-Mu_M(i+1,:))*...
                                (2*[VART(i,1),VART(i,2);VART(i,2),VART(i,3)])^-1*...
                                (Mu_D(i+1,:)-Mu_M(i+1,:))'...
                                + 1/2 * log(det([VART(i,1),VART(i,2);VART(i,2),VART(i,3)]));
                        end
                    end
                    
                    Xi2 = Xi2(Data(i_Salt).N_Cells(2:end)~=0,:);
                    Objective_Mu = sum(sum(Xi2));
                    
                    Sig_M = Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Sig;
                    
                    N_Times = length(obj.tt);
                    Objective_Sig = zeros(1,N_Times);
                    p = size(Data(i_Salt).Scatter,2);
                    for i_time = 1:N_Times
                        n = Data(i_Salt).N_Cells(i_time);
                        if n~=0;
                            Scat_D = squeeze(Data(i_Salt).Scatter(i_time,:,:));
                            ind = 1;
                            Sig = zeros(p,p);
                            for j = 1:p
                                for k = j:p
                                    Sig(j,k) = Sig_M(i_time,ind);
                                    Sig(k,j) = Sig_M(i_time,ind);
                                    ind=ind+1;
                                end
                                Sig(j,j) = max(1e-3,Sig(j,j));
                            end
                            %                 Objective_Sig(i_time) =-(log(det(Scat_D))*((n-p-1)/2)+...
                            %                     (-trace(Sig\Scat_D)/2)-...
                            %                     log(2)*((n*p/2))-...
                            %                     log(det(Sig))*(n/2) -...
                            %                     LogmGamma(p,n/2));
                            Objective_Sig(i_time) =-(log(det(Scat_D))*((n-p-1)/2)+...
                                (-trace(Sig\Scat_D)/2)-...
                                log(2)*((n*p/2))-...
                                log(det(Sig))*(n/2) -...
                                LogmGamma(p,n/2));
                            %                 Objective_Sig(i_time) =trace(Sig\Scat_D)/2+log(det(Sig))*n/2;
                        elseif n==0
                            Objective_Sig(i_time)=0;
                        end
                    end
                    %                 Objective(i_Salt) = Objective_Mu + sum(Objective_Sig);
                    %                 Objective_Sig(2:end)
                    Objective(i_Salt) = Objective_Mu + sum(Objective_Sig(3:end));
                end
            case 0  % Means, but using Model Variance
                Mu_M = Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Mu.*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                Mu_D = Data(i_Salt).Trajectories(:,1:size(Mu_M,2)).*repmat(Data(i_Salt).N_Cells,1,size(Data(i_Salt).Output_Mu,2));
                VART = max(Data(i_Salt).Variance(2:end,:),1e-6);
                Xi2 = (Mu_D(2:end,:)-Mu_M(2:end,:)).^2./(2*VART.*repmat(Data(i_Salt).N_Cells(2:end),1,size(Data(i_Salt).Output_Mu,2)));
                Xi2 = Xi2(Data(i_Salt).N_Cells(2:end)~=0,:);
                Objective_Mu = sum(sum(Xi2));

                Sig_M = Dynamics(i_Salt).Trajectories*Data(i_Salt).Output_Sig;
                                
                N_Times = length(obj.tt);
                Objective_Sig = zeros(1,N_Times);
                p = size(Data(i_Salt).Scatter,2);
                for i_time = 1:N_Times
                    n = Data(i_Salt).N_Cells(i_time);
                    if n~=0;
                        Scat_D = squeeze(Data(i_Salt).Scatter(i_time,:,:));
                        ind = 1;
                        Sig = zeros(p,p);
                        for j = 1:p
                            for k = j:p
                                Sig(j,k) = Sig_M(i_time,ind);
                                Sig(k,j) = Sig_M(i_time,ind);
                                ind=ind+1;
                            end
                            Sig(j,j) = max(1e-3,Sig(j,j));
                        end
                        %                 Objective_Sig(i_time) =-(log(det(Scat_D))*((n-p-1)/2)+...
                        %                     (-trace(Sig\Scat_D)/2)-...
                        %                     log(2)*((n*p/2))-...
                        %                     log(det(Sig))*(n/2) -...
                        %                     LogmGamma(p,n/2));
                        Objective_Sig(i_time) =-(log(det(Scat_D))*((n-p-1)/2)+...
                            (-trace(Sig\Scat_D)/2)-...
                            log(2)*((n*p/2))-...
                            log(det(Sig))*(n/2) -...
                            LogmGamma(p,n/2));
                        %                 Objective_Sig(i_time) =trace(Sig\Scat_D)/2+log(det(Sig))*n/2;
                    elseif n==0
                        Objective_Sig(i_time)=0;
                    end
                end
%                 Objective(i_Salt) = Objective_Mu + sum(Objective_Sig);
%                 Objective_Sig(2:end)
                Objective(i_Salt) = Objective_Mu + sum(Objective_Sig(3:end));
        end
    elseif strcmp(obj.Type,'Moments_4')        
        N_Times = length(obj.tt);
        for i_time = 3:N_Times
            n = Data(i_Salt).N_Cells(i_time);
            if n~=0
                switch obj.Spatial
                    case {'Spatial','Fixed'}
                        x = Data(i_Salt).Trajectories(i_time,:)';
                        FUNNAME = ['conv_mom_',num2str(obj.MGenes*2),'_',num2str(obj.Moment_Order)];
                        if ~exist(FUNNAME,'file')
                            make_mom_conv_fun(obj.MGenes*2,obj.Moment_Order,FUNNAME)
                        end
                        funhandle = str2func(FUNNAME);
                        inds = find(max(Data(i_Salt).Bred(2:end,1:obj.Num_States),[],2)==0);
                        xm_uncen = Dynamics(i_Salt).Trajectories(i_time,inds)';
                        for i=length(xm_uncen):-1:1; v{i} = xm_uncen(i); end
                        xm = funhandle(v{:});
                        L = 1/n*[xm([3,4,6,7,8])';...
                            xm([4,5,7,8,9])';...
                            xm([6,7])', xm(10) - (n-3)/(n-1)*(xm(3))^2, xm(11)-(n-3)/(n-1)*xm(3)*xm(4),                xm(12)-(n-3)/(n-1)*xm(3)*xm(5);...
                            xm([7,8])', xm(11)-(n-3)/(n-1)*xm(3)*xm(4), xm(12)-(n-2)/(n-1)*xm(4)^2+1/(n-1)*xm(3)*xm(5),xm(13)-(n-3)/(n-1)*xm(4)*xm(5);...
                            xm([8,9])', xm(12)-(n-3)/(n-1)*xm(3)*xm(5), xm(13)-(n-3)/(n-1)*xm(4)*xm(5),                xm(14) - (n-3)/(n-1)*(xm(5))^2];
                        logdetL = log(det(L));
                        x_mod = [xm_uncen(1:2);xm(3:5)];
                        Objective_Sig(i_time) = -1/2*logdetL - 1/2*(x-x_mod)'*(L^(-1))*(x-x_mod);
                    case 'NonSpatial'
                        x = Data(i_Salt).Trajectories(i_time,:)';
                        FUNNAME = ['conv_mom_',num2str(obj.MGenes),'_',num2str(obj.Moment_Order)];
                        if ~exist(FUNNAME,'file')
                            make_mom_conv_fun(obj.MGenes,obj.Moment_Order,FUNNAME)
                        end
                        funhandle = str2func(FUNNAME);
                        inds = find(max(Data(i_Salt).Bred(2:end,1:obj.Num_States),[],2)==0);
                        xm_uncen = Dynamics(i_Salt).Trajectories(i_time,inds)';
                        for i=length(xm_uncen):-1:1
                            v{i} = xm_uncen(i); 
                        end
                        xm = funhandle(v{:});
                        L = 1/n*[xm([2,3])';xm(3),xm(4) - (n-3)/(n-1)*(xm(2))^2];
                        logdetL = log(det(L));
                        x_mod = [xm_uncen(1);xm(2)];
                        Objective_Sig(i_time) = -1/2*logdetL - 1/2*(x-x_mod)'*(L^(-1))*(x-x_mod);
                        
                end
            end
        end
        Objective(i_Salt) = -sum(Objective_Sig(3:end));
        
    elseif strcmp(obj.Type,'Distributions')
%         [Dynamics,N_Space] = hog_model_FSP_class(obj);
        N_Times = length(obj.tt);
        err_ar = zeros(1,N_Times);
        Data(i_Salt).Trajectories(1,N_Space(2)+1,N_Space(3)+1)=0;
        for i_time = 1:N_Times
            if strcmp(obj.Spatial,'Spatial')||strcmp(obj.Spatial,'Fixed')
                for j_nuc=1:N_Space(2)+1
                    for j_cyt=1:N_Space(3)+1
                        err_ar(i_time) = err_ar(i_time)+...     % Likelihood
                            (Data(i_Salt).Trajectories(i_time,j_nuc,j_cyt)*...
                            log(Dynamics(i_time,j_nuc,j_cyt)));
                    end
                end
            elseif strcmp(obj.Spatial,'NonSpatial')&&PARS.Num_Genes==2
                for j_Spec1=1:N_Space(2)+1
                    for j_Spec2=1:N_Space(3)+1
                        err_ar(i_time) = err_ar(i_time)+...     % Likelihood
                            (Data(i_Salt).Trajectories(i_time,j_Spec1,j_Spec2)*...
                            log(Dynamics(i_time,j_Spec1,j_Spec2)));
                    end
                end
            elseif strcmp(obj.Spatial,'NonSpatial')
                A = squeeze(Data(i_Salt).Trajectories(i_time,:));
                B = squeeze(Dynamics(i_time,:,1));
                for i=1:length(B)
                    err_ar(i_time) = err_ar(i_time)+A(i)*log(B(i)); % Likelihood
                end
            end
        end
%         Objective(i_Salt) = -sum(err_ar);  % This converts the objective to the minus log-likelihood (a positive number that we seek to minimize).
        Objective(i_Salt) = -sum(err_ar(2:end));  % This converts the objective to the minus log-likelihood (a positive number that we seek to minimize).
        %     Objective = Objective+sum(abs(Par_Vector(abs(Par_Vector)>1e6)));  %added to restrict parameters from getting too large.
        Dynamics=[];
    end
end
Objective=sum(Objective);

    