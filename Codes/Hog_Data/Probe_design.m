close all
PD_STL1 = xlsread('Hog_Data/STL1_CTT1_probes.xlsx','STL1');
PD_STL1 = PD_STL1(:,2);

for i=1:1630
    C(i) = sum(PD_STL1-1630+i>0);
end

stairs(C); hold on;

%
PD_CTT1 = xlsread('Hog_Data/STL1_CTT1_probes.xlsx','CTT1');
PD_CTT1 = PD_CTT1(:,2);

for i=1:1630
    C(i) = sum(PD_CTT1-1630+i>0);
end

stairs(C)

save('Hog_Data/Hog_Probes.mat','PD_*')