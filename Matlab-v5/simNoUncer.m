function [ y, c, a, v ] = simNoUncer(policyA1,EV,startingA)

% This function takes the policy functions and value functions in an environment
% where there is no uncertainty, along with starting assets and returns 
% simulated paths of income, consumption assets and value

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to
global T r 
global Agrid  Ygrid numSims interpMethod 


%% ------------------------------------------------------------------------
% Initialise arrays that will hold the paths of income consumption, value
% and assets

% Arguments for output
y = NaN(T, numSims);            % income
c = NaN(T, numSims);            % consumption
v = NaN(T, numSims);            % value
a = NaN(T + 1,numSims);         % this is the path at the start of each period, so we include the 'start' of death

%% ------------------------------------------------------------------------
% Obtain paths using the initial condition and the policy and value
% functions
%-------------------------------------------------------------------------%

for s = 1:1:numSims                   % loop through individuals
    a(1, s) = startingA;   
    for t = 1:1:T                     % loop through time periods for a particular individual
        y(t  , s)   = Ygrid(t);                        
        v(t  , s)   = interp1(Agrid(t, :),EV(t, :),a(t, s),interpMethod, 'extrap');                               
        a(t+1, s)   = interp1(Agrid(t, :), policyA1(t, :) ,a(t, s),interpMethod, 'extrap');        
        c(t  , s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+r));
    end   %t      
end % s
  
 
 
end