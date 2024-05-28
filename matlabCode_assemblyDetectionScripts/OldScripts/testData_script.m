load('test_data.mat');
nneuTestData=size(spM,1); 
BinSizesTestData=[0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4 0.6 0.85 1.5];
MaxLagsTestData=[10 10 10 10 10 10 10 10 10 10 10];
[assemblyTestData]=Main_assemblies_detection(spM,MaxLagsTestData,BinSizesTestData);
[As_across_binsTestData,As_across_bins_indexTestData]=assemblies_across_bins(assemblyTestData,BinSizesTestData);
[AmatrixTestData,BinvectorTestData,Unit_orderTestData,As_orderTestData]=assembly_assignment_matrix(As_across_binsTestData, nneuTestData, BinSizesTestData, 'raw');