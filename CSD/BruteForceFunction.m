function [currentBest_Kb, currentBest_Kg, currentBest_MSE, currentBest_ModelData_storage, t_total] = BruteForceFunction(currentExperimentalData, b, g, FAILSAFE_MAXIMUM, iterations, gridSize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Matlab able to read functions in subdirectory of this folder.

current_folder = "D:\Software\MATLAB\Scripts\P&O,CIT\CSD";
addpath(current_folder + "\Experimental")
addpath(current_folder + "\Model")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in experimental data
%{
% This is the path to the folder in which there are folders within folders
% with photos.
path_data = "D:\P&O\Microscopy\Photos\";

% Define all values which should not be included in the final data. These
% can be made with the getSampleName function.
exceptions = ["3003-Sample1-4.tif"; "3003-Sample1-22.tif"];

% Define all days, these are the names of the folders of the Sample data.
% Could be done automatically, but this is more time efficient in the end.
days = ["2803", "2903", "3003"];

% Gives width of desired bin size in micrometers
desiredBinSize = 5;

% Load all experimental data in equal bins.
[equalDataStorage, newBins] = CSD_Experimental(path_data, exceptions, days, desiredBinSize);

% Find the middle of every new bin.
average_newBins = AverageEqualBins(newBins);

for k = 1: length(equalDataStorage)
    equalDataStorage{k} = equalDataStorage{k} ./ average_newBins;
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First iteration. This one is totally different from other loop, because
% at this moment in time there is no interval in which a minimum would be
% found. So this tries to find some interval in which there can be worked
% later.



% Not yet iterating over all experimental data, so I just take the first
% data set. Can be upped to 2 and 3, but I don't yet do this. This is good
% for testing if everything works.



% Initiating all things that need to be saved if there is any improvement.
% This is the best indices (row, column), which are the iter_Kb and
% iter_Kg. These are needed further down the road.
% The best Kb and Kg and the best set of data are saved as well. The best
% ModelData set is saved, because I can easily check if it is actually a
% numerically stable solution.
currentBest_iter_Kb = 1;
currentBest_Kb = 0;

currentBest_iter_Kg = 1;
currentBest_Kg = 0;

currentBest_ModelData = [];
currentBest_ModelData_storage = {};

% Keeps count how many times the currentBest values change. Needed to store
% the vector of the model.
currentBest_counter = 1;


% This is going to be to check if the Mean Squared Error is improved. This
% needs to be a ridiculuously high value, so that there is at least one
% value that will be better in the set of values tried.
currentBest_MSE = 10E50;


% DO NOT CHANGE. 
% How many iterations should be done on Kb and Kg. The model is build on 
% powers of 10, so very little unnecessary iterations are done. 
iter_Kb_max = 10;
iter_Kg_max = 10;


% These are exit conditions. Not all of them are really functional at the
% moment, but because the model works now I am not going to change them.
FAILSAFE = 0;
FLAG = true;
CHECKFORSOLVER = true;

% After this many times making a grid of iter_Kb_max by iter_Kg_max the
% loop stops. This is to prevent infinite looping.
%FAILSAFE_MAXIMUM = 10; % About a minute or maybe two normally. Depends a bit on the used computer.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting kinetic values of the problem. At the moment there is no loop
% over g and b, but this could be implemented. So at this moment the loop
% only goes over Kb and Kg.
% These are just some start values. From literature I have a vague idea
% what to expect, but the computer needs to find the specific values for
% this comparison.
% Birth and growth constants [variable]
Kb_start = 1E15;
Kg_start = 1E-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first loop

tic % start timing

% Loop over all experimental data (Preselected at the moment)
%for experimentalIndex = 1: length(equalDataStorage)
    
    % Select a current experimental dataset to work with
    %currentExperimentalData = equalDataStorage{experimentalIndex};


    % Start the loop for the first iteration. Stops if one of the following
    % conditions is met:
    % - Iteration ran too many times (FAILSAFE)
    % - The minimum is not located at the edges of the MSE_storage-matrix,
    %   thus a minimum is found. (FLAG)
    % - Numeric instability when Kg and Kb get too close to oneanother (CHECKFORSOLVER)
    while (FAILSAFE < FAILSAFE_MAXIMUM & FLAG & CHECKFORSOLVER)

          % Initiate storages. MSE_storage will be used later on. The
          % others can be removed for efficiency later on, but this is not
          % the thing that is going to slow down the code.
          MSE_storage = zeros(iter_Kb_max, iter_Kg_max);
          Kb_storage = zeros(iter_Kb_max, iter_Kg_max);
          Kg_storage = zeros(iter_Kb_max, iter_Kg_max);
        
        % Iteration over Kb
        for iter_Kb = 1: iter_Kb_max 
        
            % Birth constant [variable]
            % The minus one implies that the left boundary value is also
            % used. The right boundary isn't, but that is not a problem,
            % because that means that in iteration 1 values 1-9 are
            % evaluated and in iteration 2 10-19. So this value will be
            % found if it is the one we look for.
            % Helper function is found on the last lines of this .m file.
            Kb = UpKNotStart(Kb_start, (iter_Kb-1));
            

            % Iteration over Kg
            for iter_Kg = 1:iter_Kg_max %iteration over K_g
            
                % Growth constant [variable]
                Kg = UpKNotStart(Kg_start, (iter_Kg - 1));
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Check if loop stands still
                % Remove once the code works properly.
                %[iter_Kb, iter_Kg] 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Store current values
                Kb_storage(iter_Kb, iter_Kg) = Kb;
                Kg_storage(iter_Kb, iter_Kg) = Kg;
                
                % Calculate the model CSD for the current kinetic
                % parameters. This works exactly the same as the main.m
                % with the code in it. It is now a function of the kinetic
                % parameters.
                [z, L, C, T] = CSD_model(Kg, g, Kb, b);
        
                % The code produces a figure. I can't find where, but this at least closes it.
                close all 
                

                % Select only the last curve. This would be according to
                % the model the CSD that came out of the pipe, so this is
                % the one that should be compared to the experimental data.
                modelFinalCurve = L(end, :)*10^(-9);
                        

                % Shrink model vector with one element. At the moment it takes the
                % average of two edges. How or what it does is not of importance at
                % this moment. The code runs stuck somewhere else, so if this
                % function fills my vector with unicors I wouldn't care.
                modelData = AverageModelVector(modelFinalCurve);


                % This is the MSE between the experimental data and the
                % model data. This is manually entered, because other
                % Matlab functions don't work properly.
                model_lengths = abs(modelData);
                experimental_lengths = currentExperimentalData;

                current_MSE = CalculateMSEVectors(model_lengths, experimental_lengths);
        
                % Save Mean Squared Error
                MSE_storage(iter_Kb, iter_Kg) = current_MSE;
                

                % If the Mean Squared Error currently is better than the
                % best one I got already, means that my model now fits
                % better than the last best model. So save everything as
                % the best. 
                if current_MSE <= currentBest_MSE
                   currentBest_iter_Kb = iter_Kb;
                   currentBest_Kb = Kb;
        
                   currentBest_iter_Kg = iter_Kg;
                   currentBest_Kg = Kg;
        
                   currentBest_ModelData = modelData;
                   currentBest_ModelData_storage{currentBest_counter} = currentBest_ModelData;
                    
                   currentBest_MSE = current_MSE;
                   currentBest_counter = currentBest_counter + 1;

                end

            end % end for iteration over K_g
        end % end for iteration over K_b
            
                
        % This whole block of ifs and elses checks if the best value is
        % located against the side of the matrix. If this is the case then
        % move the matrix, because I cannot select values that are located
        % above or below that value. It needs to not lay against the side
        % of a matrix for a proper result further on. This could be
        % improved, but at the moment of writing it works as expected.
        if ~(currentBest_iter_Kb == 1) && ~(currentBest_iter_Kb == iter_Kb_max)
            if ~(currentBest_iter_Kg == 1) && ~(currentBest_iter_Kg == iter_Kg_max)
                        FLAG = false;
                        disp("A proper minimum is found in the first iteration for b " + string(b) + " and g " +string(g))
                        disp("This message was produced at: " + string(datetime('now')))
            elseif currentBest_iter_Kg == 1
                        Kg_start = Lowerk_start(Kg_start, iter_Kg_max);
            elseif currentBest_iter_Kg == iter_Kg_max
                        Kg_start = UpK_start(Kg_start, iter_Kg_max);
            end
        
        elseif currentBest_iter_Kb == 1
            Kb_start = Lowerk_start(Kb_start, iter_Kb_max);
        
            if currentBest_iter_Kg == 1
                    Kg_start = Lowerk_start(Kg_start, iter_Kg_max); 
            elseif currentBest_iter_Kg == iter_Kg_max
                    Kg_start = UpK_start(Kg_start, iter_Kg_max);
            end
        
        else
            Kb_start = UpK_start(Kb_start, iter_Kb_max);
        
            if currentBest_iter_Kg == 1
                        Kg_start = Lowerk_start(Kg_start, iter_Kg_max);
            elseif currentBest_iter_Kg == iter_Kg_max
                        Kg_start = UpK_start(Kg_start, iter_Kg_max);
            end
        end
        

        % Check for numeric stability. This problem is situated in the
        % ode15s solver and if not adressed, the solver runs infinitely.
        if Kb_start/Kg_start < 10E5
            CHECKFORSOLVER = false;
            disp("Manual input needed. Solving interval is too small. Kg_start and Kb_start are too close too one another. This results in infinite waiting time.")
            currentBest_Kb = "Does not converge";
            currentBest_Kg = "Does not converge";
            currentBest_MSE = "Does not converge";
            currentBest_ModelData_storage = "Does not converge";
            t_total = "Does not converge";
        end
        
        % One iteration has passed.
        FAILSAFE = FAILSAFE + 1;
        
        % Give feedback to the user if the code is stopped forcefully.
        if FAILSAFE == FAILSAFE_MAXIMUM
            disp("Failsafe is used. Too many iterations required.")
            currentBest_Kb = "Does not converge";
            currentBest_Kg = "Does not converge";
            currentBest_MSE = "Does not converge";
            currentBest_ModelData_storage = "Does not converge";
            t_total = "Does not converge";
        end
    
    
    end % end for while

%end


t_iterationOne = toc; % end timing
%beep

if ~(FAILSAFE == FAILSAFE_MAXIMUM) & (CHECKFORSOLVER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other iterations.
% From now on the interval around the value is taken. Divided in small
% steps an iterated over. This is done for as many times as the user wants.
% The grid is square and divides into squares. This does not have to be the
% case, but I don't see any reason why not to. It is not more efficient or
% convenient I believe.

% Specify amount of iterations
%iterations = 4;


% Specify in how many parts the interval needs to be split.
%gridSize = 20;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Take the values just beside the found minimum. The minimum is located
    % somewehere in this interval. 
    intervalKb = [UpKNotStart(Kb_start, ((currentBest_iter_Kb-1)- 1 ) ), UpKNotStart(Kb_start, ((currentBest_iter_Kb-1) + 1))];
    intervalKg = [UpKNotStart(Kg_start, ((currentBest_iter_Kg-1) - 1) ), UpKNotStart(Kg_start, ((currentBest_iter_Kg-1) + 1))];
    
    
    % Sort the values from smallest to greatest. Necessary for iteration in a
    % moment.
    intervalKb = sort(intervalKb);
    intervalKg = sort(intervalKg);
    
    
    % Find span of iteration. Necessary to decide how big the step is going to
    % be.
    intervalKb_span = abs((intervalKb(2) - intervalKb(1)));
    intervalKg_span = abs((intervalKg(2) - intervalKg(1)));
    
    
    % Define step size
    intervalKb_step = intervalKb_span/gridSize;
    intervalKg_step = intervalKg_span/gridSize;
    
    
    tic % start timing
    
    % Start iterating 
    for k = 1: iterations
        MSE_storage = zeros(iter_Kb_max, iter_Kg_max);
        
        % k
    
        % Start iterating over Kb
        for iter_Kb = 1:gridSize
    
    
            % This is the reason of sorting. The sign used is dependant on if
            % it is ascending or not. This assumes that the minimum is not an
            % edge value. This is a reasonable assumption, because otherwise
            % the minimum wouldn't have been somewhere in the middle of the
            % interval in the first iteration.
            Kb = intervalKb(1) + intervalKb_step *(iter_Kb - 1);
    
    
            % Start iterating over Kg
            for iter_Kg = 1:gridSize
                
                Kg = intervalKg(1) + intervalKg_step *(iter_Kg - 1);
    
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Only for checking. Should be removed once code is finished.
                % Just to see on which iteration we are, because this sometimes
                % takes a lot of time.
                %[k, iter_Kb, iter_Kg]
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                % Code from here until first end is exactly the same as in
                % iterations 1.
                [z, L, C, T] = CSD_model(Kg, g, Kb, b);
            
                close all 
                
                modelFinalCurve = L(end, :)*10^(-9);
                           
                modelData = AverageModelVector(modelFinalCurve);
                
    
                model_lengths = abs(modelData);
                experimental_lengths = currentExperimentalData;
    
                current_MSE = CalculateMSEVectors(model_lengths, experimental_lengths);
            
    
                MSE_storage(iter_Kb, iter_Kg) = current_MSE;
                    
    
                
                if current_MSE < currentBest_MSE
                   currentBest_iter_Kb = iter_Kb;
                   currentBest_Kb = Kb;
            
                   currentBest_iter_Kg = iter_Kg;
                   currentBest_Kg = Kg;
            
                   currentBest_ModelData = modelData;
                   currentBest_ModelData_storage{currentBest_counter} = currentBest_ModelData; 
                    
                   currentBest_MSE = current_MSE;
                   currentBest_counter = currentBest_counter + 1;
                end
    
            end % end Kg iteration
        end % end Kb iteration
        
        % Update the interval for Kb
        newLeftEdge_Kb = intervalKb(1) + intervalKb_step *((currentBest_iter_Kb - 1) - 1);
        newRightEdge_Kb = intervalKb(1) + intervalKb_step *((currentBest_iter_Kb - 1) + 1);
    
        % This should not be sorted, because the first interval was sorted an
        % the manner in which this was constructed. The sort() function is not
        % very efficient, thus omitting it will cost a second. (On 100
        % iterations that will save quite some time.)
        intervalKb = [newLeftEdge_Kb, newRightEdge_Kb];
    
        % Update the interval for Kg
        newLeftEdge_Kg = intervalKg(1) + intervalKg_step *((currentBest_iter_Kg - 1) - 1);
        newRightEdge_Kg = intervalKg(1) + intervalKg_step *((currentBest_iter_Kg - 1) + 1);
    
    
        intervalKg = [newLeftEdge_Kg, newRightEdge_Kg];
        
        % Update the span of the interval. abs() function no longer needed.
        % Increases efficiency minutely
        intervalKb_span = (intervalKb(2) - intervalKb(1));
        intervalKg_span = (intervalKg(2) - intervalKg(1));
    
        % Update the step
        intervalKb_step = intervalKb_span/gridSize;
        intervalKg_step = intervalKg_span/gridSize;
    
    end % end k (main) iteration
    
    
    t_iterationOther = toc; % Stop timer
    
    t_total = t_iterationOne + t_iterationOther;
end
%{
for i = 1:5
    beep
    pause(2)
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions

% iter is still in this function, because it was used in previous possible
% minimum-finding-algorithms. All these functions just prevent that I have
% to change 50 lines if I decide to change the algorithm. 
function output = UpK_start(startInput, iter)
    power10Input = floor(log10(2*startInput));
    output = 10^(floor(power10Input))*10 - 10^(power10Input) ;
end

function output = Lowerk_start(startInput, iter)
    power10Input = floor(log10(startInput)); 
    output = 10^(power10Input)/10 + 10^(floor(log10(startInput) - 1));
end

function output = UpKNotStart(InputStart, iter)
    output = InputStart + (10^(floor(log10(InputStart))))*(iter);
end

% Take all the values of the bins and takes the average of every bin.
% Assumes equal bins.
function averageBins = AverageEqualBins(binsInput)
    stepSize = binsInput(2) - binsInput(1);
    leftEdge = (binsInput(2) + binsInput(1)) /2;
    rightEdge = (binsInput(end) + binsInput(end-1))/2;

    averageBins = [leftEdge: stepSize: rightEdge];
end

function output = CalculateMSEVectors(inputVector1, inputVector2)
    output = sqrt( (mean(inputVector1(:) - inputVector2(:))).^2 );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end