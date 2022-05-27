function [RMSE] = CalculateRMSE(inputVector1, inputVector2)
%CalculateRMSE calculates the Root Mean Square error between to vectors and
%returns the RMSE value.

    RMSE = sqrt( (mean(inputVector1(:) - inputVector2(:))).^2 );

end