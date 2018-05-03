clear all
Ns =5;
n = 4;
make_mom_conv_fun(Ns,n,['conv_mom_',num2str(Ns),'_',num2str(n),'.m'])

%%
% t = sym('t',[Ns 1],'real');
% syms m100 m010 m001 m200 m110 m101 m020 m011 m002 m300 m210 m201 m120 m111 m102 m030 m021 m012 m003 'real'
% uncent = [m100 m010 m001 m200 m110 m020 m101 m011 m002 m300 m210 m120 m030 m201 m111 m021 m102 m012 m003]'
% for i=1:length(uncent); v{i} = uncent(i); end
% cent = bbbb(v{:})
% 
% %%
% 
% u = cent(1:Ns);
% k=Ns;
% for i=1:Ns
%     for j = i:Ns
%         k=k+1;
%         s(i,j) = cent(k);
%         s(j,i) = cent(k);
%     end
% end
% 
% mgf = exp(u'*t+1/2*t'*s*t);
% al = [1 1 0];
% M_al = subs(diff(diff(diff(mgf,al(1),'t1'),al(2),'t2'),al(3),'t3'),t,[0;0;0])
% 
% 
%%
% 
% syms m1 m2 m11 m12 m22 m111 m112 m122 m222 'real'
% uncent = [1 m1 m2 m11 m12 m22 m111 m112 m122 m222]'
% cent = uncent2cent(uncent,Ns,3)
