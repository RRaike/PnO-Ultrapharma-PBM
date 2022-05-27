function [equalDataStorage, newBins, all_days_bin_sums, bins] = CSD_Experimental(path_data, exceptions, days, desiredBinSize)
% CSD_Experimental is the function form of the Test_CSD file. This is used
% automatic comparison between the model_CSD and the CSD_Experimental. It
% returns the equally sized bins and data. No figures will be made in this
% file, because this would only slow this program.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% bins, thus could be defined outside of the loop. These bins are made by
% Fiji, thus this is not something chosen by myself.
bins = [0: 10: 190, 200: 25: 375, 400: 50: 500];
bins_average = [5: 10: 195, 212.5: 25: 387.5, 425: 50: 475];


% Loop over all the data over all days
for k = 1: length(all_days_bin_sums)

    % Take the data of a day.
    bins_sum = all_days_bin_sums{k};
    
    % The program of newData and newBins, should probably be split in
    % separate files, because the newBins is always the same in this case.
    % Yet, this program can work for all strange bin sizes thrown at it. As
    % long as desiredBinSize and all bins are divisable. 
    [newData, newBins] = EqualizeBins(bins_sum, bins, desiredBinSize);

    equalDataStorage{k} = newData;
end

end