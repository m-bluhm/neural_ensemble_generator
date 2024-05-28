load('simData1_noGroups');
nneuSim1_ng=size(neurontimematrix,1); 
BinSizes_sim1_ng=[0.015 0.025 0.04 0.06 0.085 0.15 0.25];
MaxLags_sim1_ng=[10 10 10 10 10 10 10];
[assembly_sim1_ng]=Main_assemblies_detection(neurontimematrix,MaxLags_sim1_ng,BinSizes_sim1_ng);
[As_across_bins_sim1_ng,As_across_bins_index_sim1_ng]=assemblies_across_bins(assembly_sim1_ng,BinSizes_sim1_ng);
[Amatrix_sim1_ng,Binvector_sim1_ng,Unit_order_sim1_ng,As_order_sim1_ng]=assembly_assignment_matrix(As_across_bins_sim1_ng, nneuSim1_ng, BinSizes_sim1_ng, 'raw');