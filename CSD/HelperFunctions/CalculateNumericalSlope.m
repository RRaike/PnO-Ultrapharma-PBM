function [leftToStart, rightToStart, slope] = CalculateNumericalSlope(startInput, stepSize)
% CalculateNumericalSlope gives back the value of the value to the left of
% the starting point. The value to the right of it and the average slope
% going from the left to the right point.


leftToStart = startInput - stepSize;
rightToStart = startInput + stepSize;

slope = (rightToStart - leftToStart)/(2* stepSize);

end