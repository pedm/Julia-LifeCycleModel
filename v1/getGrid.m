function [grid] = getGrid(minongrid, maxongrid, GridPoints, method)

%NB: This code borrows some code from Chris Carroll of John's Hopkins University
%His website is here http://www.econ2.jhu.edu/people/ccarroll/ and his
%notes on solving dynamic models are very good - they can be found
%here: http://www.econ2.jhu.edu/people/ccarroll/SolvingMicroDSOPs.pdf


%% ------------------------------------------------------------------------ 
% Explanation

% We need to a grid from a to b. A basic approach could involve spacing out
% the point equally. The following line of code would achieve this:
%     grid= linspace(a, b, GridPoints);

% If we want to space them out so that the growth rate of the distance
% between spaces is equal we achieve this by spacing out the logs of grid 
% equally. We could achieve this by the following line of code:
%     grid= exp(linspace(log(a), log(b), GridPoints));

% This is problematic if a<=0. The approach taken here is to get a grid
% from log(1) to log(b - a + 1), exponentiate and then subrtact 1
% so that we have a grid from 0 to b-a. If we have add a to each point we
% get the desired result - a log-spaced grid from a to b



span = maxongrid - minongrid;     % b - a                  

if strcmp(method, 'equalsteps')
    grid= linspace(minongrid, span, GridPoints);
elseif strcmp(method, 'logsteps')
  loggrid = linspace(log(1), log(1+span), GridPoints);
  grid = exp(loggrid)-1;
elseif strcmp(method, '3logsteps')
  loggrid = linspace(log(1+log(1+log(1))), log(1+log(1+log(1+span))), GridPoints);
  grid = exp(exp(exp(loggrid)-1)-1)-1;   
elseif strcmp(method, '5logsteps')
  loggrid = linspace(log(1+log(1+log(1+log(1+log(1))))), log(1+log(1+log(1+log(1+log(1+span))))), GridPoints);
  grid = exp(exp(exp(exp(exp(loggrid)-1)-1)-1)-1)-1;   
elseif strcmp(method, '10logsteps')
  loggrid = linspace(log(1+log(1+log(1+log(1+log(1+log(1+log(1+log(1+log(1+log(1)))))))))), log(1+log(1+log(1+log(1+log(1+log(1+log(1+log(1+log(1+log(1+span)))))))))), GridPoints);
  grid = exp(exp(exp(exp(exp(exp(exp(exp(exp(exp(loggrid)-1)-1)-1)-1)-1)-1)-1)-1)-1)-1;   
else
    error('Error in getgrid. You have entered an invalid method for choosing the distance between grid points. Method must be one of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps.');
end


grid = grid + minongrid*ones(1, GridPoints);

end

