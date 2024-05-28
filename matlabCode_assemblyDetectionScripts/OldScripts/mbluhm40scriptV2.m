rawpath = '\\middcloud.middlebury.edu\middfiles\SpecialProjects\DashResearch\InfraslowProject\max\SimulatedDataComplete\40secs\';
figpath = '\\middcloud.middlebury.edu\middfiles\SpecialProjects\DashResearch\InfraslowProject\max\SimulatedDataComplete\40secs\figures\';
outputpath = '\\middcloud.middlebury.edu\middfiles\SpecialProjects\DashResearch\InfraslowProject\max\SimulatedDataComplete\40secs\outputFiles\';

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
   
  % Iterate over the cells of the bin field
    % Assuming assembly is your 1x1 struct
    assembly_detected = false; % Initialize the flag to false

    % Iterate over the cells of the bin field
    all_empty = true; % Assume all cells are empty
    for j = 1:numel(assembly.bin)
        if ~isempty(assembly.bin{j}) % If any cell is not empty
            all_empty = false; % Update the flag
            break; % Exit the loop since we found a non-empty cell
        end
    end

    % If all cells are not empty, at least one assembly is detected
    if ~all_empty
        disp('At least one assembly is detected.');
        [As_across_bins, As_across_bins_index] = assemblies_across_bins(assembly, BinSizes);
        
        %display parameter for unpruned assemblies
            %'raw': displayed as returned by main_assemblies_detection, in the order of the temporal precision associated with them
            % 'ordunit': units are rearranged in order to separate assembly from non-assembly units 
            % 'clustered': assemblies and participating units are clustered such that assemblies with similar unit compositions are shown next to each other, but not necessarily in the order of their temporal precisions
        display = 'raw';
        % Visualizing - all initial assemblies found
        [Amatrix, Binvector, Unit_order, As_order] = assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display);
        out_file = fullfile(figpath, ['detected_assembly_raw_', erase(nameArray{i}, '.csv'), '.png']); % Modify output filename as needed
        fig=gcf;  
        saveas(fig,out_file)
        close all
        
        % Pruning using the specified criteria
        criteria = 'biggest'; % Change this to the desired criteria
        % criteria = 'distance';
        [pruned_assemblies, pruned_assemblies_idx] = pruning_across_bins(As_across_bins, As_across_bins_index, nneu, criteria);
        % [pruned_assemblies, pruned_assemblies_idx] = pruning_across_bins(As_across_bins, As_across_bins_index, nneu, criteria ,'pvalue');

        % Visualizing - pruned assemblies for the specified criteria
        [Amatrix, Binvector, Unit_order, As_order] = assembly_assignment_matrix(pruned_assemblies, nneu, BinSizes, display);
        out_file_pruned= fullfile(figpath, ['detected_assembly_pruned_', criteria, '_', erase(nameArray{i}, '.csv'), '.png']); % Modify output filename as needed
        fig=gcf; 
        saveas(gcf,out_file_pruned)
        close all
        
        % Computing activity for pruned assemblies - specified criteria
        [assembly_activation_full] = Assembly_activity_function(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'full');
        [assembly_activation_partial] = Assembly_activity_function(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'partial');
        [assembly_activation_combined] = Assembly_activity_function(pruned_assemblies, assembly, spM2, BinSizes, 'beginning', 'combined');
        compare_activity = [assembly_activation_full{1}(:,1), assembly_activation_full{1}(:,2), assembly_activation_partial{1}(:,2), assembly_activation_combined{1}(:,2)];
     
         % Saving output for the specified criteria
        output_filename = fullfile(outputpath, ['output_pruned_', criteria, '_', erase(nameArray{i}, '.csv')]);
        save(output_filename, 'assembly', 'assembly_activation_partial', 'assembly_activation_full', 'assembly_activation_combined', 'compare_activity', 'pruned_assemblies', 'pruned_assemblies_idx', 'spM2', 'nneu');
    else
        disp('No assembly is detected.');
        filename = ['noResult_', erase(nameArray{i}, '.csv'), '.txt'];
        filepath = fullfile(figpath, filename);
        fileID = fopen(filepath, 'w');
        fprintf(fileID, [nameArray{i}, ': No assemblies are detected.\n']);
        fclose(fileID);
        
        filepath = fullfile(outputpath, filename);
        fileID = fopen(filepath, 'w');
        fprintf(fileID, [nameArray{i}, ': No assemblies are detected.\n']);
        fclose(fileID);
    end
end
