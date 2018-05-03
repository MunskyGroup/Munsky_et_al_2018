function [v] = map_factors(x,K,pref)
v=[];
if nargin<3
    pref = [];
end
for i=1:K
    if ~isempty(x(i).x)
        if isstruct(x(i).x)
            tmp = map_factors(x(i).x,K,[pref,i]);
            for l=1:length(tmp)
                v(end+1).v = tmp(l).v;
            end
        elseif x(i).x==1
            v(end+1).v = sort([pref,i]);
        end
    end
end
