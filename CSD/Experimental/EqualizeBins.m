function [newData, newBins] = EqualizeBins(dataInput, binsInput, desiredBinSize)
% EqualizeBins distributes values in a bin equally over smaller bins.
% This assumes that the desiredBinSize and the existing bins are divisable.

if ~((length(dataInput) + 1) == length(binsInput))
    disp("Size of given data and bins are not compatible. " + ...
         "Result from EqualizeBins cannot be trusted")
end

% Find the maximal bin size in binsInput.
maxBinSize = max(binsInput);

% Split the existing span of bins in equally sized bins.
newBins = [0:desiredBinSize: maxBinSize];

% NewData is like the InputData and InputBins. The length of this matrix,
% should be one less than the matrix of the bins. This optimizes the
% program insted of letting this value grow. The bins vector is also 1 row,
% so this will be as well. 
newData = zeros(1, length(newBins) - 1);

% Index is needed for correct assignment of data in matrix.
IndexNewData = 1;

% Because the first bin is zero and thus not interesting.
% Loop over all existing bins.
for k = 1:(length(binsInput)-1)

    % Find the left edge of the current bin
    leftEdgeBin = binsInput(k);

    % Find the right edge of the current bin
    rightEdgeBin = binsInput(k+1);
    
    % Assuming that the rightEdge is always bigger than the leftEdge value.
    binWidth = rightEdgeBin - leftEdgeBin;

    % Amount of new bins needed to reach the same edges.
    amountNewBins = binWidth/desiredBinSize;
    
    % Split dataInput of set bin evenly over new bins.
    dataNewBin = dataInput(k)/amountNewBins;
    
    % Assign value to new bins.
    newData(IndexNewData: IndexNewData + amountNewBins - 1) = dataNewBin;
    
    % Update current index
    IndexNewData = IndexNewData + amountNewBins;
end

% Convert newBins to micrometers, which is its actual unit.
newBins = newBins*10^-6; 

end