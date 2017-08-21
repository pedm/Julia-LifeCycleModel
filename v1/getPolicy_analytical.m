function [ policyA1, policyC] = getPolicy_analytical

% This function gets analytical policy functions in a world where
% a) there is no uncertainty and 
% b) borrowing is allowed up to the natural borrowing constraint

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to

global T beta r gamma 
global numPointsA Agrid 

%% ------------------------------------------------------------------------ 
% Initialise output arrays
policyC     = NaN(T, numPointsA);
policyA1    = NaN(T, numPointsA);

%% ------------------------------------------------------------------------ 
% The folowing is best read with the notes

alpha = (beta^(1/gamma))*((1+r)^((1-gamma)/gamma));


for ixt = 1:1:T
    periodsLeft = T- ixt + 1;    
    for ixA = 1:1:numPointsA
        
        if (abs(alpha - 1) < 1e-5) 
            policyC(ixt, ixA) =  Agrid(ixt, ixA) / periodsLeft;
        else           
            policyC(ixt, ixA) = ((1-alpha) / (1- (alpha^periodsLeft))) * Agrid(ixt, ixA);
        end
                   
        policyA1(ixt, ixA) = (1 + r) * (Agrid(ixt, ixA) - policyC(ixt, ixA));
    end    
end
        
end

