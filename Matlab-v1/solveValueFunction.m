function [ policyA1, policyC, V ] = solveValueFunction

%This function obtains the value function for each time period and
%the policy function (i.e. optimal next-period asset choice) for each time
%period. From there we can work the optimal consumption level.

%The approach taken is by backwards recursion. The optimisation each period
%is carried out using 'fminbnd'. This is an in-built optimiser in Matlab.
%The optimisation routine it uses is known as the 'golden search method'

global T r tol minCons
global numPointsA Agrid 
global V1 Agrid1                % tomorrow's value function and asset grid

%% ------------------------------------------------------------------------ 
% GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

% Matrices to hold the policy and value functions 
V        = NaN(T+1, numPointsA);
policyA1 = NaN(T,   numPointsA);
policyC  = NaN(T,   numPointsA);        


%% ------------------------------------------------------------------------ 
%Set the terminal value function to 0
V(T + 1,:) = 0; 

%% ------------------------------------------------------------------------ 
% SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
% BACKWARDS TO ZERO, ONE PERIOD AT A TIME

for ixt=T:-1:1                                 % Loop from time T-1 to 1
    V1  = V(ixt + 1,:) ;                       %  Get tomorrow's value function
    Agrid1 = Agrid(ixt + 1, :);                % Get tomorrow's asset grid
    for ixA = 1:1:numPointsA                   % points on asset grid
        
        % Solve problem at grid points in assets      
        % ---------------------------------------------------------                        
            % Information for optimisation
            A    = Agrid(ixt, ixA);            % assets today
            lbA1 = Agrid1(1);                  % lower bound: assets tomorrow
            ubA1 = (A - minCons)*(1+r);        % upper bound: assets tomorrow
            
            % Compute solution 
            if (ubA1 - lbA1 < minCons)                         % if liquidity constrained
                negV = objectivefunc(lbA1, A); 
                policyA1(ixt,ixA) = lbA1;
            else                                                 % if interior solution
                [policyA1(ixt,ixA), negV] = ...
                    fminbnd(@(A1) objectivefunc(A1, A), lbA1, ubA1, optimset('TolX',tol));                     
            end % if (ubA1 - lbA1 < minCons)          

            % Store solution and its value
            policyC(ixt, ixA) = A - policyA1(ixt, ixA)/(1+r);
            V(ixt, ixA)       = -negV; 
    end %ixA
            
    fprintf('Passed period %d of %d.\n',ixt, T)
end %ixt

end %function

