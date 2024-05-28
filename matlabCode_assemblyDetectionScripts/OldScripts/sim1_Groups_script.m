load('simData1_Groups');
nneuSim1_g=size(neurontimematrix2,1); 
BinSizes_sim1_g=[0.015 0.025 0.04 0.06 0.085 0.15 0.25];
MaxLags_sim1_g=[5 5 5 5 5 5 5];
[assembly_sim1_g]=Main_assemblies_detection(neurontimematrix2,MaxLags_sim1_g,BinSizes_sim1_g);
[As_across_bins_sim1_g,As_across_bins_index_sim1_g]=assemblies_across_bins(assembly_sim1_g,BinSizes_sim1_g);
[Amatrix_sim1_g,Binvector_sim1_g,Unit_order_sim1_g,As_order_sim1_g]=assembly_assignment_matrix(As_across_bins_sim1_g, nneuSim1_g, BinSizes_sim1_g, 'raw');
[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins_sim1_g,As_across_bins_index_sim1_g,nneuSim1_g,'biggest');
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins_pr, nneuSim1_g, BinSizes_sim1_g, 'raw');