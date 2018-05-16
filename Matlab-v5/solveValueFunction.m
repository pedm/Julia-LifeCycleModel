function [ policyA1, policyC, V, EV, EdU ] = solveValueFunction

%This function obtains the value function for each time period and
%the policy function (i.e. optimal next-period asset choice) for each time
%period. From there we can work the optimal consumption level.

%The approach taken is by backwards recursion. The optimisation each period
%is carried out using 'fminbnd'. This is an in-built optimiser in Matlab.
%The optimisation routine it uses is known as the 'golden search method'

global T r tol minCons
global numPointsA numPointsY Agrid Ygrid incTransitionMrx 
global Agrid1 EV1

%% ------------------------------------------------------------------------ 
% GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

% Matrices to hold the policy, value and marginal utility functions 
V        = NaN(T+1, numPointsA, numPointsY);
policyA1 = NaN(T,   numPointsA, numPointsY);
policyC  = NaN(T,   numPointsA, numPointsY);        
dU       = NaN(T,   numPointsA, numPointsY);

%Matrices to hold expected value and marginal utility functions 
EV  = NaN(T+1, numPointsA, numPointsY);
EdU = NaN(T,   numPointsA, numPointsY);


%% ------------------------------------------------------------------------ 
%Set the terminal value function and expected value function to 0

EV(T + 1, :,:)  = 0;          % continuation value at T-1
V(T + 1,:,:) = 0; 
%% ------------------------------------------------------------------------ 
% SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
% BACKWARDS TO ZERO, ONE PERIOD AT A TIME

for ixt=T:-1:1                               % Loop from time T-1 to 1
    Agrid1 = Agrid(ixt + 1, :);               % The grid on assets tomorrow
    
    for ixA = 1:1:numPointsA                   % points on asset grid

        
        % STEP 1. solve problem at grid points in assets and income
        % ---------------------------------------------------------
        for ixY = 1:1:numPointsY               % points on income grid
                        
            % Value of income and information for optimisation
            A    = Agrid(ixt, ixA);            % assets today
            Y    = Ygrid(ixt, ixY);           % income today
            lbA1 = Agrid(ixt + 1, 1);          % lower bound: assets tomorrow
            ubA1 = (A + Y - minCons)*(1+r);              % upper bound: assets tomorrow
            EV1  = EV(ixt + 1,:, ixY);         % relevant section of EV matrix (in assets tomorrow)
            
            % Compute solution 
            if (ubA1 - lbA1 < minCons)                         % if liquidity constrained
                negV = objectivefunc(lbA1, A, Y); 
                policyA1(ixt,ixA,ixY) = lbA1;
            else                                             % if interior solution
                [policyA1(ixt,ixA,ixY), negV] = ...
                    fminbnd(@(A1) objectivefunc(A1, A, Y), lbA1, ubA1, optimset('TolX',tol));                                                                                                                                      
            end % if (ubA1 - lbA1 < minCons)         

            % Store solution and its value
            policyC(ixt, ixA, ixY) = A + Y - policyA1(ixt, ixA, ixY)/(1+r);
            V(ixt, ixA, ixY)       = -negV; 
            dU(ixt, ixA, ixY)      = getmargutility(policyC(ixt, ixA, ixY));
            
        end %ixY


        % STEP 2. integrate out income today conditional on income
        % yesterday to get EV and EdU
        % --------------------------------------------------------
        realisedV(:,:) = V(ixt, ixA, :);
        realiseddU(:,:) = dU(ixt, ixA, :);
        for ixY = 1:1:numPointsY
            EV(ixt, ixA, ixY)  = incTransitionMrx(ixY,:)*realisedV;
            EdU(ixt, ixA, ixY) = incTransitionMrx(ixY,:)*realiseddU;
        end %ixY
                   
    end %ixA

    fprintf('Passed period %d of %d.\n',ixt, T)
    
end %ixt


%Check the stochastic discount factors
checkSDF(policyA1, EV, dU);

end %function

