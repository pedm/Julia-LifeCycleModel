%Check various inputs

%%
% Check only one of solveUsingValueFunction and solveUsingEulerEquation is
% equal to 1
if (solveUsingValueFunction * solveUsingEulerEquation ~= 0 && solveUsingValueFunction + solveUsingEulerEquation  ~= 1)
    error('Exactly 1 of solveUsingValueFunction and solveUsingEulerFunction must be set equal to 1')
end


%%
% Check if uncertainty is off that numPointsY is set equal to 1
if (numPointsY~=1) && (isUncertainty == 0)
    error('There is no uncertainty but numPointsY was not set to 1.')
end

%%
% Check that the value of rho is not greater than 1. Warn if close to 1
if (rho>0.999) || (rho <-0.999)
    if (rho>1)  || (rho < 1)     
    error('rho is greater than 1. This code solves a stationary income process. rho greater than 1 implies a non-stationary process')
    else
    warning('rho is greater than 0.99. This code solves a stationary income process. As rho gets closer to 1 the process becomes nonstationary - possibility of numerical instability')        
    end
end

%%
% Check that the standard deviation is not too small
if (sigma<1e-10) && (isUncertainty == 1)
    if (sigma<=0) 
        error('sigma is less than or equal to zero')
    else
        warning('sigma is very small and close to zero - possibility of numerical instability. Consider turning uncertainty off')
    end
end



%% Check a number of inputs that need to be either 0 or 1
if (linearise ~= 0) && (linearise ~=1)
    error('linearise should be either 0 or 1')
end

if (borrowingAllowed ~= 0) && (borrowingAllowed ~=1)
    error('borrowingAllowed should be either 0 or 1')
end

if (isUncertainty ~= 0) && (isUncertainty ~=1)
    error('isUncertainty should be either 0 or 1')
end

