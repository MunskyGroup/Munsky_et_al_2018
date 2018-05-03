function [Data] = Process_Hog_Data_Class(obj)
%% Load Data.
load(obj.Data_File);
Data.Times = obj.tt'*60; % times in seconds.

if (strcmp(obj.Model_Obj.Pars_Type,'Vector_Given')||strcmp(obj.Model_Obj.Pars_Type,'Vector'))&&obj.MGenes==2
    Total_Gene_States = obj.Model_Obj.NStates^2;
elseif strcmp(obj.Model_Obj.Pars_Type,'Connected')||obj.MGenes==1;
    Total_Gene_States = obj.Model_Obj.NStates;
end

switch obj.Spatial
    case {'Spatial','Fixed'}
        N_spec =  Total_Gene_States+2*obj.MGenes;
    case 'NonSpatial'
        N_spec =  Total_Gene_States+obj.MGenes;
end
for i_Salt = 1:length(obj.Salt)
    switch obj.Salt(i_Salt)
        case 0.4
            IExp = 2;
        case 0.2
            IExp = 1;
    end
    
    mRNA_g1_rep1_nuc = RNAFISHexp(IExp).gene(1).rep(1).DAPIth(2).type(2).RNAraw;
    mRNA_g1_rep1_cyt = RNAFISHexp(IExp).gene(1).rep(1).DAPIth(2).type(3).RNAraw;
    mRNA_g2_rep1_nuc = RNAFISHexp(IExp).gene(2).rep(1).DAPIth(2).type(2).RNAraw;
    mRNA_g2_rep1_cyt = RNAFISHexp(IExp).gene(2).rep(1).DAPIth(2).type(3).RNAraw;
    
    if IExp==2  % Use the 1st and 3rd data sets.
        mRNA_g1_rep2_nuc  = RNAFISHexp(IExp).gene(1).rep(2).DAPIth(2).type(2).RNAraw;
        mRNA_g1_rep2_cyt  = RNAFISHexp(IExp).gene(1).rep(2).DAPIth(2).type(3).RNAraw;
        mRNA_g2_rep2_nuc  = RNAFISHexp(IExp).gene(2).rep(2).DAPIth(2).type(2).RNAraw;
        mRNA_g2_rep2_cyt  = RNAFISHexp(IExp).gene(2).rep(2).DAPIth(2).type(3).RNAraw;
        mRNA_g1_rep3_nuc(:,[1,3:16])  = RNAFISHexp(IExp).gene(1).rep(3).DAPIth(2).type(2).RNAraw;
        mRNA_g1_rep3_cyt(:,[1,3:16])  = RNAFISHexp(IExp).gene(1).rep(3).DAPIth(2).type(3).RNAraw;
        mRNA_g2_rep3_nuc(:,[1,3:16])  = RNAFISHexp(IExp).gene(2).rep(3).DAPIth(2).type(2).RNAraw;
        mRNA_g2_rep3_cyt(:,[1,3:16])  = RNAFISHexp(IExp).gene(2).rep(3).DAPIth(2).type(3).RNAraw;
        mRNA_g1_rep3_nuc(:,2) = NaN;  % In this case, we are missing data for the second time point.
        mRNA_g1_rep3_cyt(:,2) = NaN;
        mRNA_g2_rep3_nuc(:,2) = NaN;
        mRNA_g2_rep3_cyt(:,2) = NaN;
        
    else  % Use the 1st and 2nd data sets.
        mRNA_g1_rep2_nuc = RNAFISHexp(IExp).gene(1).rep(2).DAPIth(2).type(2).RNAraw;
        mRNA_g1_rep2_cyt = RNAFISHexp(IExp).gene(1).rep(2).DAPIth(2).type(3).RNAraw;
        mRNA_g2_rep2_nuc = RNAFISHexp(IExp).gene(2).rep(2).DAPIth(2).type(2).RNAraw;
        mRNA_g2_rep2_cyt = RNAFISHexp(IExp).gene(2).rep(2).DAPIth(2).type(3).RNAraw;
    end
    
    Nt = size(mRNA_g1_rep2_nuc,2);
    for iT = 1:Nt
        
        switch obj.Reps
            case '_A'
                TMP = [mRNA_g1_rep1_nuc(:,iT),mRNA_g1_rep1_cyt(:,iT),...
                    mRNA_g2_rep1_nuc(:,iT),mRNA_g2_rep1_cyt(:,iT)];
            case '_B'
                TMP = [mRNA_g1_rep2_nuc(:,iT),mRNA_g1_rep2_cyt(:,iT),...
                    mRNA_g2_rep2_nuc(:,iT),mRNA_g2_rep2_cyt(:,iT)];
            case '_C'
                if IExp~=2
%                     disp('This experiment does not have a 3rd replicate. Using (B) instead of (C)')
                    TMP = [mRNA_g1_rep2_nuc(:,iT),mRNA_g1_rep2_cyt(:,iT),...
                        mRNA_g2_rep2_nuc(:,iT),mRNA_g2_rep2_cyt(:,iT)];
                else
                    TMP = [mRNA_g1_rep3_nuc(:,iT),mRNA_g1_rep3_cyt(:,iT),...
                        mRNA_g2_rep3_nuc(:,iT),mRNA_g2_rep3_cyt(:,iT)];
                end
            case '_AB'
                TMP = [mRNA_g1_rep1_nuc(:,iT),mRNA_g1_rep1_cyt(:,iT),...
                    mRNA_g2_rep1_nuc(:,iT),mRNA_g2_rep1_cyt(:,iT);...
                    mRNA_g1_rep2_nuc(:,iT),mRNA_g1_rep2_cyt(:,iT),...
                    mRNA_g2_rep2_nuc(:,iT),mRNA_g2_rep2_cyt(:,iT)];
            case '_AC'
                if IExp~=2
%                     disp('This experiment does not have a 3rd replicate. Using (A,B) instead of (A,C).')
                    TMP = [mRNA_g1_rep1_nuc(:,iT),mRNA_g1_rep1_cyt(:,iT),...
                        mRNA_g2_rep1_nuc(:,iT),mRNA_g2_rep1_cyt(:,iT);...
                        mRNA_g1_rep2_nuc(:,iT),mRNA_g1_rep2_cyt(:,iT),...
                        mRNA_g2_rep2_nuc(:,iT),mRNA_g2_rep2_cyt(:,iT)];
                else
                    TMP = [mRNA_g1_rep1_nuc(:,iT),mRNA_g1_rep1_cyt(:,iT),...
                        mRNA_g2_rep1_nuc(:,iT),mRNA_g2_rep1_cyt(:,iT);...
                        mRNA_g1_rep3_nuc(:,iT),mRNA_g1_rep3_cyt(:,iT),...
                        mRNA_g2_rep3_nuc(:,iT),mRNA_g2_rep3_cyt(:,iT)];
                end
        end
        
        Ind_Cells = find(~isnan(sum(TMP,2)));
        
        switch obj.Spatial
            case {'Spatial','Fixed'}
                mRNA(iT).DATA = TMP(Ind_Cells,:);
            case 'NonSpatial'
                mRNA(iT).DATA = TMP(Ind_Cells,[1:2:end])+TMP(Ind_Cells,[2:2:end]);
        end
        
        switch obj.Type
            case {'Moments','Moments_Plus','Moments_Corrected','Moments_4'}
                mu(iT,:) = mean(mRNA(iT).DATA);% Now to assign the means
                
                % Here we are computiong an unbiased variance estimator
                SIGMA2 = cov(mRNA(iT).DATA,0);
                
                % Now to assign the second moments
                switch obj.Spatial
                    case {'Spatial','Fixed'}
                        sig2(iT,1:4) = SIGMA2(1,1:4);
                        sig2(iT,5:7) = SIGMA2(2,2:4);
                        sig2(iT,8:9) = SIGMA2(3,3:4);
                        sig2(iT,10) = SIGMA2(4,4);
                    case 'NonSpatial'
                        sig2(iT,1) = SIGMA2(1,1);
                        sig2(iT,2) = SIGMA2(1,2);
                        sig2(iT,3) = SIGMA2(2,2);
                end
                
                Data(i_Salt).Scatter(iT,:,:) = (mRNA(iT).DATA-repmat(mu(iT,:),size(mRNA(iT).DATA,1),1))'*...
                    (mRNA(iT).DATA-repmat(mu(iT,:),size(mRNA(iT).DATA,1),1));
                
                for i=1:size(Data(i_Salt).Scatter,2)
                    Data(i_Salt).Scatter(iT,i,i) = max(1e-3,Data(i_Salt).Scatter(iT,i,i));
                end
                
                Data(i_Salt).N_Cells(iT,1) = size(mRNA(iT).DATA,1);
                
            case {'Distributions','Distributions_TS'}
                switch obj.Spatial
                    case {'Spatial','Fixed'}
                        if strcmp(obj.Gene,'STL1')
                            Nuc_STL1_Max = max(mRNA(iT).DATA(:,1));
                            Cyt_STL1_Max = max(mRNA(iT).DATA(:,2));
                            Hists_STL1 = zeros(Nuc_STL1_Max+1,Cyt_STL1_Max+1);
                            for i = 0:Nuc_STL1_Max
                                for j = 0:Cyt_STL1_Max
                                    Hists_STL1(i+1,j+1) = sum((mRNA(iT).DATA(:,1)==i).*(mRNA(iT).DATA(:,2)==j));
                                end
                            end
                            Data(i_Salt).Trajectories(iT,1:Nuc_STL1_Max+1,1:Cyt_STL1_Max+1) = Hists_STL1;
                        elseif strcmp(obj.Gene,'CTT1')
                            Nuc_CTT1_Max = max(mRNA(iT).DATA(:,3));
                            Cyt_CTT1_Max = max(mRNA(iT).DATA(:,4));
                            Hists_CTT1 = zeros(Nuc_CTT1_Max+1,Cyt_CTT1_Max+1);
                            for i = 0:Nuc_CTT1_Max
                                for j = 0:Cyt_CTT1_Max
                                    Hists_CTT1(i+1,j+1) = sum((mRNA(iT).DATA(:,3)==i).*(mRNA(iT).DATA(:,4)==j));
                                end
                            end
                            Data(i_Salt).Trajectories(iT,1:Nuc_CTT1_Max+1,1:Cyt_CTT1_Max+1) = Hists_CTT1;
                        elseif strcmp(obj.Gene,'STL1_and_CTT1')
                            Nuc_CTT1_Max = max(mRNA(iT).DATA(:,3));
                            Cyt_CTT1_Max = max(mRNA(iT).DATA(:,4));
                            Hists_CTT1 = zeros(Nuc_CTT1_Max+1,Cyt_CTT1_Max+1);
                            for i = 0:Nuc_CTT1_Max
                                for j = 0:Cyt_CTT1_Max
                                    Hists_CTT1(i+1,j+1) = sum((mRNA(iT).DATA(:,3)==i).*(mRNA(iT).DATA(:,4)==j));
                                end
                            end
                            Data(i_Salt).CTT1Spat.Trajectories(iT,1:Nuc_CTT1_Max+1,1:Cyt_CTT1_Max+1) = Hists_CTT1;
                            
                            Nuc_STL1_Max = max(mRNA(iT).DATA(:,1));
                            Cyt_STL1_Max = max(mRNA(iT).DATA(:,2));
                            Hists_STL1 = zeros(Nuc_STL1_Max+1,Cyt_STL1_Max+1);
                            for i = 0:Nuc_STL1_Max
                                for j = 0:Cyt_STL1_Max
                                    Hists_STL1(i+1,j+1) = sum((mRNA(iT).DATA(:,1)==i).*(mRNA(iT).DATA(:,2)==j));
                                end
                            end
                            Data(i_Salt).STL1Spat.Trajectories(iT,1:Nuc_STL1_Max+1,1:Cyt_STL1_Max+1) = Hists_STL1;
                            
                            STL1_Max = max(mRNA(iT).DATA(:,1)+mRNA(iT).DATA(:,2));
                            CTT1_Max = max(mRNA(iT).DATA(:,3)+mRNA(iT).DATA(:,4));
                            Hists_STL1_CTT1 = zeros(STL1_Max+1,CTT1_Max+1);
                            for i = 0:STL1_Max
                                for j = 0:CTT1_Max
                                    Hists_STL1_CTT1(i+1,j+1) = sum(((mRNA(iT).DATA(:,1)+mRNA(iT).DATA(:,2))==i).*...
                                        ((mRNA(iT).DATA(:,3)+mRNA(iT).DATA(:,4))==j));
                                end
                            end
                            Data(i_Salt).STL1_andCTT1.Trajectories(iT,1:STL1_Max+1,1:CTT1_Max+1) = Hists_STL1_CTT1;
                        end
                    case 'NonSpatial'
                        switch obj.Gene
                            case 'STL1'
                                STL1_Max = max(mRNA(iT).DATA(:,1));
                                Hists_STL1 = zeros(STL1_Max+1,1);
                                for i = 0:STL1_Max
                                    Hists_STL1(i+1,1) = sum(mRNA(iT).DATA(:,1)==i);
                                end
                                Data(i_Salt).Trajectories(iT,1:STL1_Max+1) = Hists_STL1;
                            case 'CTT1'
                                CTT1_Max = max(mRNA(iT).DATA(:,2));
                                Hists_CTT1 = zeros(CTT1_Max+1,1);
                                for i = 0:CTT1_Max
                                    Hists_CTT1(i+1,1) = sum(mRNA(iT).DATA(:,2)==i);
                                end
                                Data(i_Salt).Trajectories(iT,1:CTT1_Max+1) = Hists_CTT1;
                            case 'STL1_and_CTT1'
                                STL1_Max = max(mRNA(iT).DATA(:,1));
                                CTT1_Max = max(mRNA(iT).DATA(:,2));
                                Hists_STL1_CTT1 = zeros(STL1_Max+1,CTT1_Max+1);
                                for i = 0:STL1_Max
                                    for j = 0:CTT1_Max
                                        Hists_STL1_CTT1(i+1,j+1) = sum((mRNA(iT).DATA(:,1)==i).*(mRNA(iT).DATA(:,2)==j));
                                    end
                                end
                                Data(i_Salt).Trajectories(iT,1:STL1_Max+1,1:CTT1_Max+1) = Hists_STL1_CTT1;
                        end
                end
        end
        Data(i_Salt).Num_cells(iT) = size(mRNA(iT).DATA,1);
    end
    
    
    if strcmp(obj.Type(1:7),'Moments')
        switch obj.Spatial
            case {'Spatial','Fixed'}
                if obj.MGenes==1;
                    switch obj.Gene
                        case 'STL1'
                            if obj.Moment_Order==1
                                Data(i_Salt).Trajectories = mu(:,[1,2]);
                                Data(i_Salt).Variance = (sig2(:,[1,5]));
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                                Data(i_Salt).Output_Mu = zeros(N_spec,2);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+2,2) = 1;
                                Data(i_Salt).Mom_Types = {'m','m'};
                            elseif obj.Moment_Order>=2
                                Data(i_Salt).Trajectories = [mu(:,1:2),sig2(:,1),sig2(:,2),sig2(:,5)];
                                Data(i_Salt).Variance = (sig2(:,[1,5]));
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                                Data(i_Salt).Output_Mu = zeros(N_spec+sum(1:N_spec),2);
                                Data(i_Salt).Output_Sig = zeros(N_spec+sum(1:N_spec),3);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+2,2) = 1;
                                Data(i_Salt).Output_Sig(end-2,1) = 1;
                                Data(i_Salt).Output_Sig(end-1,2) = 1;
                                Data(i_Salt).Output_Sig(end,3) = 1;
                                Data(i_Salt).Mom_Types = {'m','m','v','v','v'};
                                Data(i_Salt).Scatter = Data(i_Salt).Scatter(:,1:2,1:2);
                            end
                        case 'CTT1'
                            if obj.Moment_Order==1
                                Data(i_Salt).Trajectories = mu(:,[3,4]);
                                Data(i_Salt).Variance = (sig2(:,[8,10]));
                                Data(i_Salt).Output_Mu = zeros(N_spec,2);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+2,2) = 1;
                                Data(i_Salt).Mom_Types = {'m','m'};
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                            elseif obj.Moment_Order>=2
                                Data(i_Salt).Trajectories = [mu(:,3:4),sig2(:,8),sig2(:,9),sig2(:,10)];
                                Data(i_Salt).Variance = (sig2(:,[8,10]));
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                                Data(i_Salt).Output_Mu = zeros(N_spec+sum(1:N_spec),2);
                                Data(i_Salt).Output_Sig = zeros(N_spec+sum(1:N_spec),3);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+2,2) = 1;
                                Data(i_Salt).Output_Sig(end-2,1) = 1;
                                Data(i_Salt).Output_Sig(end-1,2) = 1;
                                Data(i_Salt).Output_Sig(end,3) = 1;
                                Data(i_Salt).Mom_Types = {'m','m','v','v','v'};
                                Data(i_Salt).Scatter = Data(i_Salt).Scatter(:,3:4,3:4);
                            end
                    end
                elseif obj.MGenes==2;
                    if obj.Moment_Order==1
                        Data(i_Salt).Trajectories = mu(:,1:4);
                        Data(i_Salt).Variance = (sig2(:,[1,5,8,10]));
                        Data(i_Salt).Output_Mu = zeros(N_spec,4);
                        Data(i_Salt).Output_Mu(Total_Gene_States+1,1) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+2,2) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+3,3) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+4,4) = 1;
                        Data(i_Salt).Mom_Types = {'m','m','m','m'};
                        Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                    elseif obj.Moment_Order>=2
                        Data(i_Salt).Trajectories = [mu(:,1:4),sig2(:,1:10)];
                        Data(i_Salt).Variance = (sig2(:,[1,5,8,10]));
                        Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                        Data(i_Salt).Output_Mu = zeros(N_spec+sum(1:N_spec),4);
                        Data(i_Salt).Output_Sig = zeros(N_spec+sum(1:N_spec),10);
                        Data(i_Salt).Output_Mu(Total_Gene_States+1,1) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+2,2) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+3,3) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+4,4) = 1;
                        for i = 9:-1:0
                            Data(i_Salt).Output_Sig(N_spec+sum(1:N_spec)-i,10-i) = 1;
                        end
                        Data(i_Salt).Mom_Types = {'m','m','m','m','v','v','v','v','v','v','v','v','v','v'};
                    end
                end
            case {'NonSpatial'}
                if obj.MGenes==1;
                    switch obj.Gene
                        case 'STL1'
                            if obj.Moment_Order==1
                                Data(i_Salt).Trajectories = mu(:,1);
                                Data(i_Salt).Variance = (sig2(:,1));
                                Data(i_Salt).Output_Mu = zeros(N_spec,1);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Mom_Types = {'m'};
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                            elseif obj.Moment_Order>=2
                                Data(i_Salt).Trajectories = [mu(:,1),sig2(:,1)];
                                Data(i_Salt).Variance = (sig2(:,[1]));
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                                Data(i_Salt).Output_Mu = zeros(N_spec+sum(1:N_spec),1);
                                Data(i_Salt).Output_Sig = zeros(N_spec+sum(1:N_spec),1);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Output_Sig(end,1) = 1;
                                Data(i_Salt).Mom_Types = {'m','v'};
                                Data(i_Salt).Scatter = Data(i_Salt).Scatter(:,1,1);
                            end
                        case 'CTT1'
                            if obj.Moment_Order==1
                                Data(i_Salt).Trajectories = mu(:,2);
                                Data(i_Salt).Variance = (sig2(:,3));
                                Data(i_Salt).Output_Mu = zeros(N_spec,1);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Mom_Types = {'m'};
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                            elseif obj.Moment_Order>=2
                                Data(i_Salt).Trajectories = [mu(:,2),sig2(:,3)];
                                Data(i_Salt).Variance = (sig2(:,[3]));
                                Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                                Data(i_Salt).Output_Mu = zeros(N_spec+sum(1:N_spec),1);
                                Data(i_Salt).Output_Sig = zeros(N_spec+sum(1:N_spec),1);
                                Data(i_Salt).Output_Mu(obj.Model_Obj.NStates+1,1) = 1;
                                Data(i_Salt).Output_Sig(end,1) = 1;
                                Data(i_Salt).Mom_Types = {'m','v'};
                                Data(i_Salt).Scatter = Data(i_Salt).Scatter(:,2,2);
                            end
                    end
                    
                    
                elseif obj.MGenes==2;
                    if obj.Moment_Order==1
                        Data(i_Salt).Trajectories = mu(:,1:2);
                        Data(i_Salt).Variance = (sig2(:,[1,3]));
                        Data(i_Salt).Output_Mu = zeros(N_spec,2);
                        Data(i_Salt).Output_Mu(Total_Gene_States+1,1) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+2,2) = 1;
                        Data(i_Salt).Mom_Types = {'m','m'};
                        Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                    elseif obj.Moment_Order>=2
                        Data(i_Salt).Trajectories = [mu(:,1:2),sig2(:,1:3)];
                        Data(i_Salt).Variance = (sig2(:,[1,3]));
                        Data(i_Salt).Variance = max(Data(i_Salt).Variance,1e-1);
                        Data(i_Salt).Output_Mu = zeros(N_spec+sum(1:N_spec),2);
                        Data(i_Salt).Output_Sig = zeros(N_spec+sum(1:N_spec),3);
                        Data(i_Salt).Output_Mu(Total_Gene_States+1,1) = 1;
                        Data(i_Salt).Output_Mu(Total_Gene_States+2,2) = 1;
                        Data(i_Salt).Output_Sig(end-2,1) = 1;
                        Data(i_Salt).Output_Sig(end-1,2) = 1;
                        Data(i_Salt).Output_Sig(end,3) = 1;
                        Data(i_Salt).Mom_Types = {'m','m','v','v','v'};
                    end
                end
        end
    end
    Data(i_Salt).Times = obj.tt*60;
end



if obj.Moment_Order==4
    n=4;
    %% Find all possible combinations of up to order n for each species.
    BPos = zeros((n+1)^N_spec,N_spec);
    B = zeros(1,N_spec);
    BPos(1,:) = B;
    for k = 1:(n+1)^N_spec-1
        B(1) = B(1)+1;
        j=1;
        while B(j)==n+1
            B(j)=0;
            B(j+1) = B(j+1)+1;
            j=j+1;
        end
        BPos(k+1,:) = B;
    end
    
    %% Reorder all possible combinations by order up to total order n.
    K = sum(sum(BPos,2)<=n);
    Bred = zeros(K,N_spec);
    Kmom=0;
    nM = ones(n,1);
    for i=0:n
        J = find(sum(BPos,2)==i);
        nM(i+1) = length(J);
        Bred(Kmom+1:Kmom+nM(i+1),:)=BPos(J,:);
        Kmom=Kmom+nM(i+1);
    end
    for i_Salt = 1:length(obj.Salt)
        Data(i_Salt).Bred = Bred;
        Data(i_Salt).nM = nM;
    end
end
    
