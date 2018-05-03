function Solutions = solve_hog_class(obj)
salts = obj.Salt;
for i=1:length(obj.Salt)
    obj.Salt = salts(i);
    switch obj.Type
        case 'Moments'
            Solutions(i).moments = hog_model_moments_class(obj);
            Solutions(i).likelihood = Compute_Likelihood(obj,Solutions(i).moments,0);
            Solutions(i).likelihoodB = Compute_Likelihood(obj,Solutions(i).moments,1);
            Solutions(i).likelihoodC = Compute_Likelihood(obj,Solutions(i).moments,2);
        case 'Moments_4'
            Solutions(i).moments = hog_model_moments_class(obj);
            Solutions(i).likelihood = Compute_Likelihood(obj,Solutions(i).moments,4);
        case 'Distributions'
            [Solutions(i).distributions,N_Space] = hog_model_FSP_class(obj);
            Solutions(i).likelihood = Compute_Likelihood(obj,Solutions(i).distributions,[],N_Space);
        case 'Distributions_TS'
            Solutions(i).distributions = hog_model_FSP_TS_class(obj);
%             Solutions(i).likelihood = Compute_Likelihood(obj,Solutions(i).distributions);
        case 'SSA'
            [SSA,TS] = hog_model_SSA_class(obj);
            Solutions(i).SSA = SSA;
            Solutions(i).TS = TS;
    end
end
