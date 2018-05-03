function [SSADAT,TS] = hog_model_SSA_class(obj)
% This function will use the parameter set (PARS) anf the number of states,
% and the mutant type (Arp8_Gcn5), salt level, and number of genes, and
% then run the program to compute the first and second moments as a
% function of time.
%% The following is used to make and plot SSA runs using the identified model.
PARS = obj.Model_Obj.PARS;
PARS.Moment_Order = 2;
[~,ModStuff] = hog_model_moments_class(obj);
[SSADAT,TS] = Get_SSADAT(ModStuff,obj.tt*60,obj.SSA_Runs,obj.gene_length,obj.k_elong);