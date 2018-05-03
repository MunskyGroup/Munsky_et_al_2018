function [x] = find_factors(Bred,b,K)
x=[];
for i=2:K
    if min(b-Bred(i,:))<0
        x(i).x=[];
    elseif max(abs(b-Bred(i,:)))==0
        x(i).x=1;
    else
        b_left = b-Bred(i,:);
        x(i).x = find_factors(Bred,b_left,K);
    end   
end


