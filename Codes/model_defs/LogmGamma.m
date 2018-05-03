function log_mG = LogmGamma(p,x)
% This function is necessary for the computation of the Wishart
% distribution for the likelihood of the measured covariance matrix.
v = x+(1-[1:p])/2;
log_mG = log(pi)*p*(p-1)/4+sum(gammaln(v));