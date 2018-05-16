function [ Z, EVbetweeenZ ] = getNormDev(mu, sigma_inc, trunc, N )
    
% This function returns two vectors:
% a) Z - a vector of (N+1) normal deviates (numbers) that divide a normal
% distribution with standard deviation sigma into N segments, 
% each of which are equiprobable
% b) EVbetweeenZ - the expected value of the random variables between those
% two points


%% -------------------------------------------------------------------------   
% Initialise the output and working arrays
    
   % Output
   Z = NaN(N + 1, 1);
   EVbetweeenZ = NaN(N, 1);


%% -------------------------------------------------------------------------   
%Finding the points that divide the standard normal into N segments 
%The first and last of these should, if we are using an actual normal
%distribution would be minus and plus infinity.  In reality we use a truncated 
%normal distribution and use -trunc and +trunc (where trunc is an argument)
    Z(1) = -trunc * sigma_inc;
    Z(N + 1) = trunc * sigma_inc;

%Now recursively get the rest of the points
    for ixi = 2:1:N
        Z(ixi) = sigma_inc * stdnorminv_manual((ixi-1)/N) + mu;
    end
    
    stdZ = (Z - mu .* ones(N + 1, 1)) ./ sigma_inc;
    
%Finding the expected value within each interval (see Adda & Cooper page
%58)
    PDF = stdnormpdf_manual(stdZ);
    
    for ixi = 1:1:N
        % this comes from approx page 56 of Adda and Cooper
        EVbetweeenZ(ixi) = N .* sigma_inc .* (PDF(ixi) - PDF(ixi + 1)) + mu;
    end 

end

