function [ value ] = objectivefunc(A1, A0)

%%
%-------------------------------------------------------------------------------%
% This function returns the following quantity:
% - (u(c) +  b V( A1))
% where c is calculated from today's assets and tomorrow's assets
%-------------------------------------------------------------------------------%

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to
global beta r interpMethod            % structural model parameters
global Agrid1 V1                        % tomorrow's asset grid and tomorrow's  value function

%% ------------------------------------------------------------------------ 
%Get tomorrow's consumption (cons), the value of left over assets (VA1) and
%total value (u(c) + b * VA1

cons = A0  - (A1)/(1+r);
VA1 = interp1(Agrid1, V1 , A1, interpMethod, 'extrap');
value = utility(cons) + beta * VA1;

%% ------------------------------------------------------------------------ 
%The optimisation routine that we will use searches for the minimum of the
%function. We want the maximum. So we multiply out function here by -1 so
%that the optimiser will fill the minimum of the negative of our function,
%i.e. the maximum of our functino 

value = - value;

end

