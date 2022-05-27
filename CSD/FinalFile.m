% This file brute forces all values. In the BruteForceFunction the file
% path might need to be changed on different devices.

clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Matlab able to read functions in subdirectory of this folder.

current_folder = "D:\Software\MATLAB\Scripts\P&O,CIT\CSD";
addpath(current_folder + "\Experimental")
addpath(current_folder + "\Model")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_data = "D:\P&O\Microscopy\Photos-Old\";
exceptions = ["3003-Sample1-4.tif"; "3003-Sample1-22.tif"];
days = ["2803"];
    %, "2903", "3003"];

[experimentalData, newBins] = LoadExperimentalData(path_data, days, exceptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FAILSAFE_MAXIMUM = 20;
iterations = 2;
gridSize = 10;

bInterval = [1:1:10];
gInterval = [3:1:10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KbBestAllExperiments = cell(length(experimentalData), 1);
KbBestCurrentExperiment = zeros(length(bInterval), length(gInterval));

KgBestAllExperiments = cell(length(experimentalData), 1);
KgBestCurrentExperiment = zeros(length(bInterval), length(gInterval));

MSEBestAllExperiments = cell(length(experimentalData), 1);
MSEBestCurrentExperiment = zeros(length(bInterval), length(gInterval));

modelDataBestAllExperiments = cell(length(experimentalData), 1);
modelDataCurrentExperiment = cell(length(bInterval), length(gInterval));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for experiment = 1: length(experimentalData)
    currentExperimentalData = experimentalData{experiment};
    
    tic
    for bIndex = 1: length(bInterval)

        b = bInterval(bIndex);

        parfor gIndex = 1: length(gInterval)
            
            g = gInterval(gIndex);

            disp("Experiment: " + string(experiment))
            
            [currentBest_Kb, currentBest_Kg, currentBest_MSE, currentBest_ModelData_storage, t_total] = ...
                    BruteForceFunction(currentExperimentalData, b, g, FAILSAFE_MAXIMUM, iterations, gridSize);
            
            disp("Loop these b and g is done. b = " + string(b) + " and g = " + string(g))
            disp("This message was produced at: " + string(datetime('now')))
            disp(" ")

            KbBestCurrentExperiment(bIndex, gIndex) = currentBest_Kb;
            KgBestCurrentExperiment(bIndex, gIndex) = currentBest_Kg;
            MSEBestCurrentExperiment(bIndex, gIndex) = currentBest_MSE;
            modelDataCurrentExperiment{bIndex, gIndex} = currentBest_ModelData_storage; 
        end
    end
    t_experiment = toc;

    disp("An experiment is done. This took: " + string(t_experiment) + ".")
    KbBestAllExperiments{experiment} = KbBestCurrentExperiment;
    KgBestAllExperiments{experiment} = KgBestCurrentExperiment;
    MSEBestAllExperiments{experiment} = MSEBestCurrentExperiment;
    modelDataBestAllExperiments{experiment} = modelDataCurrentExperiment;
end


disp("All data is analyzed. This message was made on: " + string(datetime('now')))
