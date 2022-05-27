function [spanStorage, spanSizesStorage] = DataCharacteristics(Data, newBins, percentileLeft, percentileRight)
%{
Input:
    Data: Should be a cell array with in each element the matrix with
    frequency of a certain bin.
    newBins: bins of the histogram. These should be equidistant.
    percentileLeft: Defines the percentile that should be the left edge of
    the span calculation. Unitless.
    percentileRight: Defines the percentile that should be the right edge
    of the span calculation. Unitless.

Output:
    spanStorage: Cell array with each cell containing the span of the Data
        set at the same index. The unit of this output is m_{crystal}.
    spanSizesStorage: Cell array with each cell containing a matrix. This
        matrix has 3 elements. The left percentile size, the mean size and
        the right percentile size. All units are m_{crystal}.
%}


spanSizesStorage = cell(length(Data),1);
spanStorage = cell(length(Data),1);

for dataCounter = 1: length(Data)
    currentExperimentalData = Data{dataCounter};
    
    frequency_sum = sum(currentExperimentalData);
    
    normalizedHistogram = currentExperimentalData /frequency_sum;
    
    mean = 0.5;
    
    spanPercentiles = [percentileLeft; mean; percentileRight];
    
   
    
    spanIndexStorage = zeros(length(spanPercentiles), 1);
    spanFractionStorage = zeros(length(spanPercentiles), 1);
    spanSizes = zeros(length(spanPercentiles), 1);
    
    for l = 1: length(spanPercentiles)
        currentPercentile = spanPercentiles(l);
         percentileSum = 0.00;

        k = 1;
        Flag = true;
        while k <= length(normalizedHistogram) & Flag
            currentProbability = normalizedHistogram(k);
            
            if (percentileSum + currentProbability) >= currentPercentile
                percentileIndex = k;
                fractionPercentile = (currentPercentile - percentileSum)/currentProbability;
                Flag = false;
            else
                percentileSum = percentileSum + currentProbability;
                k = k+ 1;
            end
        end
    
        spanIndexStorage(l) = percentileIndex;
        spanFractionStorage(l) = fractionPercentile;
        
    
    end
    
    
    binStep = newBins(2) - newBins(1);
    
    for k = 1: length(spanIndexStorage)
        currentIndex = spanIndexStorage(k);
        currentFraction = spanFractionStorage(k);
        
        spanSizes(k) = newBins(currentIndex) + binStep* currentFraction;
        
    end
    
    span = (spanSizes(3) - spanSizes(1))/spanSizes(2);
    
    spanSizesStorage{dataCounter} = spanSizes;
    spanStorage{dataCounter} = span;
end
end