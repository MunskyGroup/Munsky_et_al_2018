
function ERR = get_Fit_errs_FSP_TS(speed,inacl,reps,gene,TS_ALL,N,DAT,spat,type)
Hog.MOD = Hog_Model;
Hog.MOD.k_elong = speed;
switch gene
    case 'STL1'
%         Hog.MOD.gene_length = 1709;
        Hog.MOD.gene_length = 1630;
    case 'CTT1'
%         Hog.MOD.gene_length = 1688;
        Hog.MOD.gene_length = 1618;
end
switch type
    case 'Distributions'
        Hog.MOD.TS_Lim = 100;
    case {'Moments_Corrected','Means'}
        Hog.MOD.TS_Lim = 30000;
end
        
Hog.MOD.Gene = gene;
Hog.MOD.Model = '4State_k21only_noloop';
Hog.MOD.MGenes = 1;
Nmax = Hog.MOD.TS_Lim;

Nt = length(Hog.MOD.tt);% Number of time points for tS.
ERR = 0;

load('Part_Intens.mat','Partial_Intens_Map');
C = ones(1,21);
Ci = ones(1,10);
for j=3:Nmax
    C = sparse(blkdiag(C,Ci));
end
Partial_Intens_Map_Red = C*Partial_Intens_Map(1:Nmax*10+1,1:Nmax+1);

for i_nacl = inacl
    %% Step 2b -- Find the TS spots intensities (MODEL)
    switch i_nacl
        case 1
            Hog.MOD.Salt = 0.2;
        case 2
            Hog.MOD.Salt = 0.4;
    end
    
    Hog.Fit_Type = type;
    Hog.MOD.Spatial=spat;    
    DIR = ['Fit_Results/',Hog.MOD.Gene,'_0p2and0p4/',Hog.MOD.Spatial,'/',Hog.Fit_Type,'/',Hog.MOD.Model];    
    Hog.MOD.Parameter_File = [DIR,'/Best_Pars_AC.mat'];
    
    Hog.MOD.Type='Distributions_TS';
    Hog.Distributions_TS = Hog_Model.solve_hog(Hog.MOD);
    inten_dist = Partial_Intens_Map_Red*Hog.Distributions_TS.distributions(1:Nmax+1,:);
    
    for irep = reps
        for it=1:16
            if ~isempty(DAT(i_nacl,irep,it).A)
                er_ar = log(inten_dist(1:99,it)).*DAT(i_nacl,irep,it).A';
                ERR = ERR - sum(er_ar(DAT(i_nacl,irep,it).A>0));
            end
        end
        
    end
end
