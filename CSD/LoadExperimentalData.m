function [ExperimentalData, newBins] = LoadExperimentalData(path_data, days, exceptions)


current_folder = "D:\Software\MATLAB\Scripts\P&O,CIT\CSD";
addpath(current_folder + "\Experimental")

% Load in experimental data

% This is the path to the folder in which there are folders within folders
% with photos.
%path_data = "D:\P&O\Microscopy\Photos\";

% Define all values which should not be included in the final data. These
% can be made with the getSampleName function.
% exceptions = ["3003-Sample1-4.tif"; "3003-Sample1-22.tif"];

% Define all days, these are the names of the folders of the Sample data.
% Could be done automatically, but this is more time efficient in the end.
%days = ["2803", "2903", "3003"];

% Gives width of desired bin size in micrometers
desiredBinSize = 5;

% Load all experimental data in equal bins.
[equalDataStorage, newBins] = CSD_Experimental(path_data, exceptions, days, desiredBinSize);

% Find the middle of every new bin.
average_newBins = AverageEqualBins(newBins);

for k = 1: length(equalDataStorage)
    ExperimentalData{k} = equalDataStorage{k} ./ average_newBins;
end

end

function averageBins = AverageEqualBins(binsInput)
    stepSize = binsInput(2) - binsInput(1);
    leftEdge = (binsInput(2) + binsInput(1)) /2;
    rightEdge = (binsInput(end) + binsInput(end-1))/2;

    averageBins = [leftEdge: stepSize: rightEdge];
end