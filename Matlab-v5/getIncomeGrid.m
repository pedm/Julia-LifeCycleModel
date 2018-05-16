function [Ygrid, Q, minInc, maxInc] = getIncomeGrid

% A function that returns:
% 1. an income grid
% 2. A Markovian transition matrix (Q) over income realisations
% 3. A vector of minimum incomes in each year
% 4. A vector of maximum incomes in each year

%% ------------------------------------------------------------------------ 
%  Declare the global variables that will be necessary in this function
global T Tretire sigma mu rho isUncertainty        % Structural economic paramters
global normBnd numPointsY                        % Numerical methods

%% ------------------------------------------------------------------------ 


%----------------------------------------%
% Scenario where there is no uncertainty
%----------------------------------------%
if isUncertainty == 0
    y = exp(mu);                 % income is set equal to the exp of the log mean
    minInc = y;                 
    maxInc = y;
    Q = [1];                    % The transition matrix Q is simply a constant 1
                                % with prob 1 each period income is 1

%----------------------------------------%
% Scenario where there is uncertainty - income draws are log normally distributed 
%----------------------------------------%
                                
elseif isUncertainty == 1


        % First get the standard deviation of income (from sigma and rho)
        sig_inc = sigma/((1-rho^2)^0.5);

        % Split the entire normal distribution into numPointsY sections that
        % are equiprobable. The output lNormDev gives the (numPointsY + 1) 
        % points that bound the sections, the output ly gives the 
        % (numPointsY) expected value in each section
        [ lNormDev, ly ] = getNormDev(mu, sig_inc, normBnd, numPointsY );    

        %---------------------%
        %Get transition matrix Q(i, j). The prob of income j in t+1
        %conditional on income i in t
        %---------------------%
        Q = NaN(numPointsY, numPointsY);             %initialise the transition matrix
        for i = 1:1:numPointsY
            for j = 1:1:numPointsY
                hiDraw = lNormDev(j+1) - (1-rho)*mu - rho * ly(i); %highest innovation that will give us income j tomorrow
                loDraw = lNormDev(j)   - (1-rho)*mu - rho * ly(i); %lowest  innovation that will give us income j tomorrow
                Q(i,j) = stdnormcdf_manual(hiDraw/sigma) - ...
                         stdnormcdf_manual(loDraw/sigma);
            end %j
            
            %Each of the rows of Q should add up to 1. But
            %due to the truncation of the normal distribution they will
            %not. So we divide through by the sum of elements in the row
            Q(i, :) = Q(i, :) ./ sum(Q(i, :));                
        end %i

        y = exp(ly);                      % Get y from log y
        minInc = exp(-normBnd * sig_inc); %Get the minimum income in each year 
        maxInc = exp(normBnd * sig_inc);  %Get the maximum income in each year 
        
        if (y(1) < 1e-4) || (y(numPointsY) > 1e5)
            warning('Combination of sigma and rho give a very high income variance. Numerical instability possible')
        end
        
end  % if isUncertainty == 0


%----------------------------------------%
%Now get a matrix, T * numPointsY that holds the grid for each income in
%each year. Do likewise with minimum and maximum income
%----------------------------------------%
Ygrid = repmat(y', T, 1);
minInc = repmat(minInc', T, 1);          
maxInc = repmat(maxInc', T, 1);

%----------------------------------------%
%Replace these arays with zeros for all years after retirement
%----------------------------------------%
if Tretire == 0         % no work (retired at birth)
   Ygrid(:, :) = 0;
   minInc(:, :) = 0;
   maxInc(:, :) = 0;
elseif (Tretire > 0) && (Tretire <=T)  %retire at some age
    Ygrid(Tretire:T, :) = 0;
    minInc(Tretire:T, :) = 0;
    maxInc(Tretire:T, :) = 0;
end



end
