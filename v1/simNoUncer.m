function [ c, a, v ] = simNoUncer(policyA1,V,startingA)

% This function takes the policy functions and value functions in an environment
% where there is no uncertainty, along with starting assets and returns 
% simulated paths of consumption assets and value

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to
global T r 
global Agrid  interpMethod 


%% ------------------------------------------------------------------------
% Initialise arrays that will hold the paths of income consumption, value
% and assets

% Arguments for output
c = NaN(T, 1);            % consumption
v = NaN(T, 1);            % value
a = NaN(T + 1,1);         % this is the path at the start of each period, so we include the 'start' of death

%% ------------------------------------------------------------------------
% Obtain paths using the initial condition and the policy and value
% functions
%-------------------------------------------------------------------------%

a(1, 1) = startingA;   
for t = 1:1:T                     % loop through time periods for a particular individual
    v(t  , 1)   = interp1(Agrid(t, :),V(t, :),a(t, 1),interpMethod, 'extrap');                               
    a(t+1, 1)   = interp1(Agrid(t, :), policyA1(t, :) ,a(t, 1),interpMethod, 'extrap');        
    c(t  , 1) = a(t, 1)  - (a(t+1, 1)/(1+r));
end   %t      
  
 
 
end