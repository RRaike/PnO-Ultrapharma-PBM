function [sum_all_bins] = FilterExtendedSummary(input, exceptions)
% Input needs to be the Extended Summary table from the Excel file. This is
% gained via the InputExcel-function.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:length(input)

    table = input{k};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get all labels and find indices of all values that need to be excluded.
    
    % Initiation
    index_exceptions = [];
    
    % Find index for the column Slice and make all the values in it strings
    % instead of chars.
    index_column_Slice = find(strcmpi(table.Properties.VariableNames,'Slice'));
    Slice_column = string( table{:, index_column_Slice} );
    
    % Find the index for every exception in the exceptions input. This will
    % be used in the next section.
    for j = 1: length(exceptions)
        current_exception = exceptions(j);
    
        current_index = find(strcmp(Slice_column, current_exception));
        index_exceptions = [index_exceptions; current_index];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get all bins of sizes, make this into usable matrix and sum all values.
    
    % Find from when the bins start and then take the matrix from that until
    % the end. All inputs are a similar format, thus need not be changed with
    % different ExtendedSummary files.
    
    index_first_bin = find(strcmpi(table.Properties.VariableNames,'x0_10'));
    bin_data = table{:, index_first_bin:end};
    
    % Initiate empty vector, wherein summated data can be put.
    sum_bins = zeros( 1, length(bin_data(1, :)));
    
    
    % Loop to take all exceptions out of the data.
    for j = 1:length(bin_data(:,1))
        current_bin = bin_data(j,:);
       
    
        if ~ ismember(j, index_exceptions)
            sum_bins = sum_bins + current_bin;
        end
    
    end

    % Store data.
    sum_bins_storage{k} = sum_bins;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum_all_bins = zeros(1, length(sum_bins_storage{1}));

% Sum of all data of one sheet. Agglomerates data.
for k = 1 : length(sum_bins_storage)
    sum_all_bins = sum_all_bins + sum_bins_storage{k};
end




end