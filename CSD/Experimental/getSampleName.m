function [SampleName] = getSampleName(date, sample_number, number)
%getSampleName needs 3 inputs to recreate the name of a sample:
%   date: format of date should be DAYMONTH as one piece.
%   sample_number: number of the sample wherein a value is searched.
%   number: given the specific tracked number within a sample.
% This function does NOT check if this value actually exists.
% The returned name is given as a string.


SampleName = int2str(date) + "-Sample" + int2str(sample_number) + "-" + ...
             int2str(number) + ".tif";

end