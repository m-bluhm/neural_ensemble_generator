clear all

%NOTE: THE FIGPATH AND RAWPATH NEED TO BE CHANGED ONCE THE ISSUE WHERE
%MATLAB IS UNABLE TO WRITE INTO THE SimulatedDataSets FOLDER IS FIXED

% CHANGEABLE PARAMETERS
dataSeries = 640; %Which dataseries (time durations) we are working with

%change to whatever folder the data is stored in.
%to change from seeded to unseeded data, change the "singleDatasetsSeeded"
%to "singleDatasetsUnseeded" and the 'seeded%dSecondData' to
%'unseeded%dSecondData'
rawpath = sprintf('\\\\middcloud.middlebury.edu\\middfiles\\SpecialProjects\\DashResearch\\InfraslowProject\\MaxBluhmSpring2024\\SimulatedDataSets\\singleDatasetsSeeded\\seeded%dSecondData', dataSeries);

display = 'raw';
 %display parameters for unpruned assemblies
        %'raw': displayed as returned by main_assemblies_detection, in the order of the temporal precision associated with them
        % 'ordunit': units are rearranged in order to separate assembly from non-assembly units 
        % 'clustered': assemblies and participating units are clustered such that assemblies with similar unit compositions are shown next to each other, but not necessarily in the order of their temporal precisions

% Pruning using the specified criteria
criteria = 'biggest'; % Change this to the desired criteria. other criteria include 'distance'

% We have to use a smaller array of bin sizes if the time duration is short. 
% Realisiticly, we should be testing different binsizes and maxlags.
% However, I never got to that point.
% BinSizes: The vector of bin widths to be tested:
% MaxLags: For each bin size specify the maximum lag lmax (in numbers of bins) to be tested (all lags
% within [−lmax, lmax ] will be tested):
if dataSeries < 160
    BinSizes = [0.015 0.025 0.04 0.06 0.085 0.15 0.25];
    MaxLags = [10 10 10 10 10 10 10];
else
    BinSizes = [0.015 0.025 0.04 0.06 0.085 0.15 0.25 0.4];
    MaxLags = [10 10 10 10 10 10 10 10];
end   

% CHANGEABLE PARAMETERS

% figpath = sprintf('%s\\figures', rawpath); %The figures will be stored in the "figures" subfolder of rawpath
% outputpath = sprintf('%s\\output', rawpath); %The matlab output will be stored in the "output" subfolder of rawpath

figpath = 'C:\MATLAB\seededResults\figures';
outputpath = 'C:\MATLAB\seededResults\output';

% Get list of files in the directory
fileList = dir(rawpath);

% Exclude '.' and '..' entries from the list
fileList = fileList(~ismember({fileList.name}, {'.', '..'}) & ~[fileList.isdir]);

% Extract filenames from the fileList structure
nameArray = {fileList.name};

for i = 1:length(nameArray)
    nameArray(i)
    clearvars -except BinSizes MaxLags rawpath nameArray figpath outputpath dataSeries criteria display i
    % Display the current file being processed

    % Get the filename
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
        %To visualize the detected assembly patterns, it is necessary to first create a new structure As_across_bins, in which assembly units are re-orderd according to their order of activation within the assembly:
        [As_across_bins, As_across_bins_index] = assemblies_across_bins(assembly, BinSizes);

        % Visualizing - all initial assemblies found
        [Amatrix, Binvector, Unit_order, As_order] = assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display);
        out_file = fullfile(figpath, ['detected_assembly_' display, '_', erase(nameArray{i}, '.csv'), '.png']); % Modify output filename as needed
        fig=gcf;  
        saveas(fig,out_file)
        close all
        %pruning
        [pruned_assemblies, pruned_assemblies_idx] = pruning_across_bins(As_across_bins, As_across_bins_index, nneu, criteria);

        % Visualizing - pruned assemblies for the specified criteria
        [Amatrix, Binvector, Unit_order, As_order] = assembly_assignment_matrix(pruned_assemblies, nneu, BinSizes, display);
        out_file_pruned= fullfile(figpath, ['detected_assembly_pruned_', criteria, '_', erase(nameArray{i}, '.csv'), '.png']); % Modify output filename as needed
        fig=gcf; 
        saveas(gcf,out_file_pruned)
        close all

        % Computing activity for pruned assemblies - specified criteria. This
        % looks at the 'full,' 'partial,' and 'combined' activation criteria
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
        try
            fileID = fopen(filepath, 'w');
            fprintf(fileID, [nameArray{i}, ': No assemblies are detected.\n']);
            fclose(fileID);
        catch
            disp('Error occurred while writing to file.');
        end

        filepath = fullfile(outputpath, filename);
        try
            fileID = fopen(filepath, 'w');
            fprintf(fileID, [nameArray{i}, ': No assemblies are detected.\n']);
            fclose(fileID);
        catch
            disp('Error occurred while writing to file.');
        end
    end
end
    
    

