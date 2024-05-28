load('neurontimematrix160secsr30groups.mat');
spM=neurontimematrix160secsr30groups;
nneu=size(spM,1); 
BinSizes=[0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4 0.6 0.85];
MaxLags=[10 10 10 10 10 10 10 10 10 10];
[assembly]=Main_assemblies_detection(spM,MaxLags,BinSizes);
[As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly,BinSizes);
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display);