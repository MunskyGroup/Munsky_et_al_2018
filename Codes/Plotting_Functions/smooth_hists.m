function [hists] = smooth_hists(hists_in,N)
%% [hists] = smooth_hists(hists_in,N)
% This function smooths the function hists_in by taking the average of the
% entries from [i-N : i+N].  
Nmax = length(hists_in);
for i=1:Nmax
    hists(i) = mean(hists_in(max(1,i-N):min(Nmax,i+N)));
end