function hog_sig = Define_Hog_Signal(t,Hog_Parameters)
% This subroputine defines the hog signal at the current instant in time.
t =max(0,t-Hog_Parameters.t_offset);
hog_sig = (1-exp(-Hog_Parameters.r1*(t)))*exp(-Hog_Parameters.r2*(t));
hog_sig = (hog_sig/(1+hog_sig*Hog_Parameters.M))^Hog_Parameters.eta/1.0727183838607e-10;