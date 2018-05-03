function [X_Array,TS] = Get_SSADAT(Model_Properties,TArray,Nruns,gene_length,k_elong)

Nx = size(Model_Properties.Stoichiometry,1); % get number of species.
NTimes = length(TArray); % get number of times.
X_Array = zeros(Nruns,Nx,NTimes);  % Preallocate trajectory.
tau_ts = linspace(0,gene_length/k_elong,1000);
parfor jrun = 1:Nruns
% for jrun = 1:Nruns
    [X_Array(jrun,:,:),TS(jrun,:,:)] = run_SSA(Model_Properties,TArray,tau_ts);
end

function [X_Array,TS] = run_SSA(Model_Properties,TArray,tau_ts)
S = Model_Properties.Stoichiometry;
Nx = size(S,1); % get number of species.
S = [zeros(Nx,1),S];
TRXNs = Model_Properties.transcription_inds+1;  % Indices of the transcription events.

W1 =  full([zeros(1,Nx);Model_Properties.W1]);
W1t = full([zeros(1,Nx);Model_Properties.Signal.W1t]);
Trate = 1;
W0 = full([Trate;Model_Properties.W0]);
W0t = full([0;Model_Properties.Signal.W0t]);
Signal = Model_Properties.Signal;

NTimes = length(TArray); % get number of times.
x0 = Model_Properties.Moments_0(1:Nx);
tstop = TArray(NTimes);


T1 = zeros(length(TArray),length(tau_ts));
for i=1:length(tau_ts)
    T1(:,i) = TArray-tau_ts(i);
end
% T1 is a matrix of earlier times. There are 1000 earlier times considered
% beginning at the time exactly one full elongation time in the past. Each
% row corresponds to one snapshot time, and each column corresponds to a
% time prior to that snapshot time.

x = x0;  %initial condition.
t = 0;   %initial time.
X_Array = zeros(length(x0),length(TArray));
TS = zeros(length(TRXNs),size(T1,1),size(T1,2));
iTime_Count = 1;
while t<tstop
    %% Time varying rates
    if strcmp(Signal.Type,'Interpolate')==1
        it =1;
        while t>Signal.Time(it)&&it<length(Signal.Time)
            it=it+1;
        end
        if it==1
            Input_Signal = Signal.Value(1);
        elseif it<=length(Signal.Time)
            Input_Signal = Signal.Value(it-1)+...
                (Signal.Value(it)-Signal.Value(it-1))*...
                (t-Signal.Time(it-1))/(Signal.Time(it)-Signal.Time(it-1));
        else
            Input_Signal = Signal.Value(length(Mean_ERK));
        end
    elseif strcmp(Signal.Type,'Function')==1
        Input_Signal = Signal.Function(t);
    end
    
    w0 = W0 + W0t*Input_Signal;   % Time varying values of W0
    w1 = W1 + W1t*Input_Signal;   % Time varying values of W1
    w = max(0,w1*x+w0); 
    %         if min(W0)<0||min(min(W1))<0
    %             w0 = sparse(max(0,W0));  %Propensities must always be positive for any State.
    %             w1 = sparse(max(0,W1));  %Propensities must always be positive for any State.
    %             SAT = norm(full(W0-sparse(max(0,W0 + W0t*1e10)))) + norm(full(W1-sparse(max(0,W1 + W1t*1e10))));
    %             if SAT~=0
    %                 w = max(0,W1*x+W0); w(1) = Trate;
    %             else
    %                 W0=W0
    %             end
    %         end
    %%
    %         w = max(0,W1*x+W0); w(1) = Trate;
    w0 = sum(w); %% Compute the sum of the prop. functions
    t = t+1/w0*log(1/rand); %% Update time of next reaction
    if t<=tstop
        while t>TArray(iTime_Count)
            X_Array(:,iTime_Count) = x;
            iTime_Count = iTime_Count+1;
        end
        r2w0=rand*w0; %% generate second random number and multiply by prop. sum
        i=1; %% initialize reaction counter
        while sum(w(1:i))<r2w0 % increment counter until sum(w(1:i)) exceeds r2w0
            i=i+1;
        end
        x = x+S(:,i); % update the configuration
        [istrxn,jtrxn]=ismember(i,TRXNs);
        if istrxn&&max((t>T1(:,end))'.*(t<TArray))~=0
            % The logic t>T1(:,end) asks if the initiation occurs recently
            % enough that that the mRNA would not have completed. The logic
            % t<tArray asks if the initiation occurs before the current
            % smFISH snapshot.
            
            TS(jtrxn,(t>T1)&(t<repmat(TArray',1,length(tau_ts)))) =...
            TS(jtrxn,(t>T1)&(t<repmat(TArray',1,length(tau_ts))))+1;
%             TS(jtrxn,(t<T1)&(t>repmat(T1(:,end),1,length(tau_ts)))&(t<repmat(TArray',1,length(tau_ts)))) =...
%             TS(jtrxn,(t<T1)&(t>repmat(T1(:,end),1,length(tau_ts)))&(t<repmat(TArray',1,length(tau_ts))))+1;
            
%         sparse(squeeze(TS))
%         [t/60 (t+max(tau_ts))/60]
%         
%         pause
%             
%             
        end
    end
end
X_Array(:,iTime_Count:end) = repmat(x,1,NTimes-(iTime_Count-1));
