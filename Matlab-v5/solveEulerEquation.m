function [ policyA1, policyC, V, EV, EdU  ] = solveEulerEquation

%% ------------------------------------------------------------------------ 
%This function obtains the value function for each time period and
%the policy function (i.e. optimal next-period asset choice) for each time period.

%The approach taken is by backwards recursion. The optimisation each period
%is carried out using 'fzero'. This function finds the zero of an arbitrary
%function. We use it to find the value of saving (and therefore
%consumption) that ensures that the Euler equation holds.

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to

global T r tol linearise minCons
global numPointsA Agrid Ygrid
global incTransitionMrx numPointsY
global Agrid1 Edu1 linEdU1 EV1      %Ass grid, expec marg util, linearised expec marg util and value tomorrow

%% ------------------------------------------------------------------------ 
% GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

% Matrices to hold the policy, value and marginal utility functions 
V        = NaN(T + 1,numPointsA, numPointsY) ;
policyA1 = NaN(T    ,numPointsA,numPointsY) ;
policyC  = NaN(T,    numPointsA, numPointsY) ;   
dU       = NaN(T,    numPointsA, numPointsY);

%Matrices to hold expected value and marginal utility functions 
EV       = NaN(T + 1,numPointsA, numPointsY);
EdU      = NaN(T,    numPointsA, numPointsY);

%A vector that will tomorrow's period's 'linearised' marginal utility 
lindU1   = NaN(numPointsA, 1);


%% ------------------------------------------------------------------------ 
%Set the terminal value function and expected value function to 0
V(T + 1,:,:)  = 0 ;
EV(T + 1,:,:)  = 0;

%% ------------------------------------------------------------------------ 
% SOLVE THE CONSUMER'S PROBLEM AT TIME T, WHEN SOLUTION IS KNOWN
% Optimal consumption is equal to assets held as there is no value to keep them for
% after death

for ixA = 1:1:numPointsA                                         % points on asset grid

    % STEP 1. solve problem for each grid point in assets and income today
    % ---------------------------------------------------------------------
    for ixY = 1:1:numPointsY                                        % points on income grid

         % Value of state variables
         Y    = Ygrid(T, ixY);                                     % income today
         A    = Agrid(T, ixA);                                      % assets today

         % Compute and store solution and its value
         policyC(T, ixA, ixY) = A + Y;                            % optimal consumption
         policyA1(T, ixA, ixY) = 0;                               % optimal next period assets
         V(T, ixA, ixY)  = utility(policyC(T, ixA, ixY));         % value of policyC
         dU(T, ixA, ixY) = getmargutility(policyC(T, ixA, ixY));  % marginal value of policyC

    end %ixY

    % STEP 2. integrate out income today conditional on income yesterday to
    % get EV and EdU 
    % ---------------------------------------------------------------------
    realisedV(:,:) = V(T, ixA, :);
    realiseddU(:,:) = dU(T, ixA, :);    
    for ixY = 1:1:numPointsY                                  % for each point on the income grid (last income)
        EV(T, ixA, ixY)  = incTransitionMrx(ixY,:)*realisedV;          % continuation value at T-1
        EdU(T, ixA, ixY) = incTransitionMrx(ixY,:)*realiseddU;         % expect marginal utility at T
    end %ixY
    
end %ixA

fprintf('Passed period %d of %d.\n',T, T)

%% ------------------------------------------------------------------------ 
% SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
% BACKWARDS TO ZERO, ONE PERIOD AT A TIME

for ixt=(T-1):-1:1      %Loop from time T-1 to 1
    Agrid1 = Agrid(ixt + 1, :);               % The grid on assets tomorrow
    for ixA = 1:1:numPointsA    %Loop through points on the asset grid

    % STEP 1. solve problem for each grid point in assets and income today
    % ---------------------------------------------------------------------

        for ixY = 1:1:numPointsY
            %   Value of income and information for optimisation
            A = Agrid(ixt, ixA);                %assets today
            Y = Ygrid(ixt, ixY);               %income today
            lbA1 = Agrid(ixt + 1, 1);           %lower bound: assets tomorrow
            ubA1 = (A + Y - minCons)*(1+r);     %upper bound: assets tomorrow
            bndForSol = [lbA1, ubA1];           %if the Euler equation has a soluton it will be within these bounds
            Edu1 = EdU(ixt + 1,:, ixY);         %relevant section of Edu matrix (in assets tomorrow)     
            if linearise == 1                   %get 'linearised' marginal utility tomorrow
                  linEdU1 = getinversemargutility(Edu1);
            end

            %   Compute solution 
            signoflowerbound = sign(eulerforzero(A, lbA1, Y));
            if (signoflowerbound == 1) || (ubA1 - lbA1 < minCons)    %if liquidity constrained 
                 policyA1(ixt, ixA, ixY) = lbA1;
            else                                                %if interior solution                              
                signofupperbound = sign(eulerforzero(A, ubA1, Y));
                if (signoflowerbound*signofupperbound == 1)
                    error('Sign of lower bound and upperbound are the same - no solution to Euler equation. Bug likely')
                end
                [policyA1(ixt, ixA, ixY)] = fzero(@(A1) eulerforzero(A, A1, Y), bndForSol,optimset('TolX',tol));
            end                                     %if (signoflowerbound == 1) 

            %   Store solution and its value          
            policyC(ixt, ixA, ixY) = A + Y - policyA1(ixt, ixA, ixY)/(1+r);
            dU(ixt, ixA, ixY) = getmargutility(policyC(ixt, ixA, ixY));
            EV1 = EV(ixt + 1,:, ixY);                %Select out the value function at each point on the asset grid tomorrow                                  
            V(ixt, ixA, ixY) = -objectivefunc(policyA1(ixt, ixA, ixY), A, Y);
        end %ixY

        % STEP 2. integrate out income today conditional on income
        % yesterday to get EV and EdU
        % --------------------------------------------------------
        realisedvalues(:,:) =V(ixt, ixA, :);        
        realisedMargUtility(:,:) = dU(ixt, ixA, :);

        for ixY = 1:1:numPointsY       
            EV(ixt, ixA, ixY) =  incTransitionMrx(ixY,:)*realisedvalues;
            EdU(ixt, ixA, ixY) =  incTransitionMrx(ixY,:)*realisedMargUtility;
        end %ixY
        
    end %ixA

    fprintf('Passed period %d of %d.\n',ixt, T)
end %ixt

end

