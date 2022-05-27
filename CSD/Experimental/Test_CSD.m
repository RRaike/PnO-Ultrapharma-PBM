clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is the path to the folder in which there are folders within folders
% with photos.
path_data = "D:\P&O\Microscopy\Photos\";

% Define all values which should not be included in the final data. These
% can be made with the getSampleName function.
exceptions = ["3003-Sample1-4.tif"; "3003-Sample1-22.tif"];

% Define all days, these are the names of the folders of the Sample data.
% Could be done automatically, but this is more time efficient in the end.
days = ["2803", "2903", "3003"];

desiredBinSize = 5;

[equalDataStorage, newBins, all_days_bin_sums, bins] = CSD_Experimental(path_data, exceptions, days, desiredBinSize);

%{
% Iteration of every day and thus over every map. There will be 1
% dataset in the form of a Excelworksheets (with multilple sheets) per day.
for k = 1: length(days)
    current_day = days(k);

    % Get the ExtendedSummary Excel file read into Matlab
    ExtendedSummary = InputExcel(path_data, current_day);

    % Read the ExtendedSummary file. Summate all the same bins in different
    % sheets. This returns the distribution of all the data of 1 day.
    % Automatically removes the exceptions of the dataset.
    current_bins_sum = FilterExtendedSummary(ExtendedSummary, exceptions);
    
    % Store the data of 1 day.
    all_days_bin_sums{k} = current_bins_sum;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define bins as defined in Excel datasheet. This are always the same
% bins, thus could be defined outside of the loop.
bins = [0: 10: 190, 200: 25: 375, 400: 50: 500];
bins_average = [5: 10: 195, 212.5: 25: 387.5, 425: 50: 475];


bins_sum = all_days_bin_sums{1};
[newData, newBins] = EqualizeBins(bins_sum, bins, 5);
%}

%{
% Initiate empty storages to store Rsquares, Means and STDs later.
Rsquare = zeros(length(days), 1);
Means = zeros(length(days), 1);
STDs = zeros(length(days), 1);

% Loop over all the data over all days
for k = 1: length(all_days_bin_sums)

    % Take the data of a day.
    bins_sum = all_days_bin_sums{k};
    
    
    % Fit distribution (Gaussian) trough the current experimental data.
    [bell, goodness] = fit(bins_average', bins_sum','gauss1');
    
    % Get the necessary data from the fit.
    current_Rsquare = goodness.rsquare;
    current_Mean = bell.b1;
    current_STD = bell.c1/sqrt(2);
    
    % Save these values. Could have been done in the last line, but this
    % increases readability I believe. Especially ones figures are being
    % made.
    Rsquare(k) = current_Rsquare;
    Means(k) = current_Mean;
    STDs(k) = current_STD;


    % Create figure
    f = figure;
    
    % Create histogram
    h = histogram();

    % Give the numerical values of every amount in a bin. 
    h.BinCounts = bins_sum;

    % Defines the bins. The bins are not equally sized, thus Matlab only
    % allows it to be draw like this.
    h.BinEdges = bins;

    % Possible Normalization of the data. Not necessary at this moment
    %h.Normalization = 'pdf';

    hold on
    
    % Draws a dot in every center of the top of a bar of the histogram.
    % Only needed for visual purposes.
    scatter(bins_average, bins_sum)
    hold on
    
    % Plot the distribution curve that is related with the data on top of
    % the histogram.
    plot(bell);

    % Define where to put to textbox. Place textbox with all relevant
    % information.
    NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
    textOnPlot = {"R^2: " + string(current_Rsquare), ...
                  "Mean: " + string(current_Mean), ...
                  "Standard deviation: " + string(current_STD)};
    annotation('textbox', [0.55 .30 .5 .6] ,'String', textOnPlot, 'FitBoxToText', 'on')
    legend('hide')
    titleName = extractBefore(days(k), 3) + "/" + extractAfter(days(k), 2) + "/2022" + ...
        " All data with a fitted distribution";
    title(titleName)
    

    % Make the name with which the figure can be saved as .png.
    NAME = days(k) + "-" + "Figure" + int2str(k) + ".png";

    % Save the file
    %saveas(f, NAME)

    
    %Kolmogorov-Smirnovtoets
end
%}

