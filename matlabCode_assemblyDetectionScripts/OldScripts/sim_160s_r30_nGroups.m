%load('neurontimematrix160secsr30noGroups');
nneu50=size(neurontimematrix160secsr30noGroups,1); 
BinSizes=[0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4 0.6 0.85];
MaxLags=[10 10 10 10 10 10 10 10 10 10];
[assembly_160s_r30_nG]=Main_assemblies_detection(neurontimematrix160secsr30noGroups,MaxLags,BinSizes);
[As_across_bins_160s_r30_nG,As_across_bins_index_160s_r30_nG]=assemblies_across_bins(assembly_160s_r30_nG,BinSizes);
[Amatrix_160s_r30_nG,Binvector_160s_r30_nG,Unit_order_160s_r30_nG,As_order_160s_r30_nG]=assembly_assignment_matrix(As_across_bins_160s_r30_nG, nneu50, BinSizes, 'raw');