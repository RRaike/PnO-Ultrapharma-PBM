function [AverageModelVector] = AverageModelVector(ModelVector)

AverageModelVector = zeros(1, length(ModelVector) - 1);

for k = 1:length(AverageModelVector)
    firstElement = ModelVector(k);
    secondElement = ModelVector(k+1);

    AverageModelVector(k) = (firstElement + secondElement)/ 2;
end

end