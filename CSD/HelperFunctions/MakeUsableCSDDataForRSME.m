function [modelData] = MakeUsableCSDDataForRSME(Kg, g, Kb, b, volume)

% Calculate model
[z, L, C, T] = CSD_model(Kg, g, Kb, b);

close all 
% Take only final curve. This would be the predicted curve at the end of
% the tube. Times the volume to go to the unit #/ m_crystal
modelFinalCurve = L(end, :) *volume;    

% Shrink the vector by one element.
modelData = AverageModelVector(modelFinalCurve);


end