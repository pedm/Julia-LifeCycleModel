function [ policyA1, policyC] = getPolicy_analytical

% This function gets analytical policy functions in a world where
% a) there is no uncertainty and 
% b) borrowing is allowed up to the natural borrowing constraint

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to

global T beta r gamma minCons 
global numPointsA Agrid Ygrid 

%% ------------------------------------------------------------------------ 
% Initialise output arrays
policyC     = NaN(T, numPointsA);
policyA1    = NaN(T, numPointsA);

%% ------------------------------------------------------------------------ 
% The folowing is best read with the notes

alpha = (beta^(1/gamma))*((1+r)^((1-gamma)/gamma));


for ixt = 1:1:T
    periodsLeft = T- ixt + 1;
    indexVec = (0:(periodsLeft-1))';                            % a vector [ 0 1       2           ... (periodsLeft-1) ]
    RtoThePowerOft = ((1+r) .* ones(periodsLeft,1)).^indexVec;  % a vector [ 1 (1+r)  (1+r)^2      ... (1 + r)^(periodsLeft-1) ]
    inverseOfRtoThePowerOft = 1./RtoThePowerOft;                % a vector [ 1 1/(1+r) 1/((1+r)^2) ... 1/((1 + r)^(periodsLeft-1)) ]

    inc = Ygrid(ixt:T, 1);                                     % income from today until death
    discountedInc = inc' * inverseOfRtoThePowerOft;             % the sum of income stream discounted by interest rate
    
    for ixA = 1:1:numPointsA
        W =  Agrid(ixt, ixA) + discountedInc;
        
        if (abs(alpha - 1) < 1e-5) 
            policyC(ixt, ixA) =  W / periodsLeft;
        else           
            policyC(ixt, ixA) = ((1-alpha) / (1- (alpha^periodsLeft))) * W;
        end
                   
        policyA1(ixt, ixA) = (1 + r) * (Agrid(ixt, ixA) + Ygrid(ixt)  - policyC(ixt, ixA));
    end    
end
        
end

