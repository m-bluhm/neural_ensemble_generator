rawpath = 'C:\Users\mbluhm\Documents\Max_Neuron_Project\seeded40SecondData';
figpath = '\\middcloud.middlebury.edu\middfiles\SpecialProjects\DashResearch\InfraslowProject\max\SimulatedDataComplete\40secsSeeded\figures';
outputpath = '\\middcloud.middlebury.edu\middfiles\SpecialProjects\DashResearch\InfraslowProject\max\SimulatedDataComplete\40secsSeeded\output';
    % nameArray = {'neuron_time_matrix_40secs_r=30_groups.csv', 'neuron_time_matrix_40secs_r=30_noGroups.csv', 'neuron_time_matrix_40secs_r=30_noSine_groups.csv', 'neuron_time_matrix_40secs_r=30_noSine_noGroups.csv', 'neuron_time_matrix_40secs_r=60_groups.csv', 'neuron_time_matrix_40secs_r=60_noGroups.csv', 'neuron_time_matrix_40secs_r=60_noSine_groups.csv', 'neuron_time_matrix_40secs_r=60_noSine_noGroups.csv'};
nameArray = {'seeded_neuron_time_matrix_40secs_r=60_noSine_groups.csv', 'seeded_neuron_time_matrix_40secs_r=60_noSine_noGroups.csv'}
BinSizes = [0.015 0.025 0.04 0.06 0.085 0.15 0.25];
MaxLags = [10 10 10 10 10 10 10];


nameArray = {'neuron_time_matrix_40secs_r=30_groups.csv', 'neuron_time_matrix_40secs_r=30_noGroups.csv', 'neuron_time_matrix_40secs_r=30_noSine_groups.csv', 'neuron_time_matrix_40secs_r=30_noSine_noGroups.csv', 'neuron_time_matrix_40secs_r=60_groups.csv', 'neuron_time_matrix_40secs_r=60_noGroups.csv', 'neuron_time_matrix_40secs_r=60_noSine_groups.csv', 'neuron_time_matrix_40secs_r=60_noSine_noGroups.csv'};
%nameArray = {'neuron_time_matrix_40secs_r=30_noGroups.csv', 'neuron_time_matrix_40secs_r=30_noSine_groups.csv', 'neuron_time_matrix_40secs_r=30_noSine_noGroups.csv', 'neuron_time_matrix_40secs_r=60_groups.csv', 'neuron_time_matrix_40secs_r=60_noGroups.csv', 'neuron_time_matrix_40secs_r=60_noSine_groups.csv', 'neuron_time_matrix_40secs_r=60_noSine_noGroups.csv'};

BinSizes = [0.015 0.025 0.04 0.06 0.085 0.15 0.25];
MaxLags = [10 10 10 10 10 10 10];

for i = 1:length(nameArray)
    nameArray(i)
    %clearvars -except BinSizes MaxLags rawpath nameArray
    filename = fullfile(rawpath, nameArray{i});
    [a, b, c] = xlsread(filename);
    spM2 = a;
    nneu = size(spM2, 1);
    
    % Assembly detection
    assembly = Main_assemblies_detection(spM2, MaxLags, BinSizes);
   
  
    % Assembly reordering
    [As_across_bins, As_across_bins_index] = assemblies_across_bins(assembly, BinSizes);
    
    % Visualizing - all initial assemblies found
    [Amatrix, Binvector, Unit_order, As_order] = assembly_assignment_matrix(As_across_bins, nneu, BinSizes, 'raw');
    
    out_file = fullfile(figpath, ['detected_assembly_raw', num2str(i),'.png']); % Modify output filename as needed
    fig=gcf;  
    saveas(fig,out_file)
    close all
    
    % Pruning
    criteria = 'biggest';
    [pruned_assemblies, pruned_assemblies_idx] = pruning_across_bins(As_across_bins, As_across_bins_index, nneu, criteria);
    
    % Visualizing - pruned assemblies
    [Amatrix, Binvector, Unit_order, As_order] = assembly_assignment_matrix(pruned_assemblies, nneu, BinSizes, 'raw');
    
    out_file = fullfile(figpath, ['detected_assembly_pruned', num2str(i),'.png']); % Modify output filename as needed
    fig=gcf; 
    saveas(gcf,out_file)
    close all
    
    % Activation using original code
    %assembly_activation_full = Assembly_activity_function(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'full');
    %assembly_activation_partial = Assembly_activity_function(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'partial');
    
    [assembly_activation_partial]=Assembly_activity_function_md(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'partial');  %now should account for multiple activations, column 2 is the average activation, column 3 is the number of ensemble activations (count) in that bin


    
    %assembly_activation_combined = Assembly_activity_function(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'combined');
    %assembly_activation_proportion = Assembly_activity_function_ebw(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'proportion');
    
    %compare_activity = [assembly_activation_full(:,1), assembly_activation_full(:,2), assembly_activation_partial(:,2), assembly_activation_combined(:,2), assembly_activation_proportion(:,2)];
    
    % Save output using the name of the CSV file
    output_filename = fullfile(rawpath, ['output_', nameArray{i}]); % Modify output filename as needed
    save(output_filename, 'assembly_activation_partial');
end





%2)  