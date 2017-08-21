%% ------------------------------------------------------------------------ 
% Dynamic Economics in Practice. 
% Monica Costa Dias and Cormac O'Dea
% Institute for Fiscal Studies
% 17th-18th September 2014

%% ------------------------------------------------------------------------ 
% DESCRIPTION
% This program solves and simulates a finite period consumption and saving 
% (cake-eating problem) problem. There is no income or uncertainty. The
% consumer starts with an endowmnent of assets and chooses how to allocate
% consumption over time. Solution is by value function iteration

%% ------------------------------------------------------------------------ 
% PREAMBLE
% Ensure that all storage spaces variables, matrices, memory, globals are 
% clean from information stored in past runs

clear all;  % clear memory
close all;  % close any graphs


%% ------------------------------------------------------------------------ 
% DECLARE VARIABLES AND MATRICES THAT WILL BE 'GLOBAL'
% explicitly set variables and matrices to be shared throughout the routine
% as globals

global beta gamma r T                               % structural model parameters 
global numPointsA Agrid                             % assets grid and dimension
global interpMethod                                 % numerical methods to be used
global tol minCons                                  % numerical constants
global plotNumber                                   % for graphing

%% ------------------------------------------------------------------------ 
% NUMERICAL METHODS
% select solution, interpolation and integration methods

interpMethod = 'linear';      % interpolation method - choose from 'linear', 'nearest', 'spline', 'pchip'

%% ------------------------------------------------------------------------ 
% NUMERICAL CONSTANTS
% set constants needed in numerical solution and simulation

% precision parameters
%--------------------------------------%
tol = 1e-10;                 % max allowed error
minCons = 1e-5;              % min allowed consumption


%% ------------------------------------------------------------------------ 
% THE ECONOMIC ENVIRONMENT
% Set values of structural economic parameters

T = 80;                      % Number of time period
r = 0.01;                   % Interest rate
beta = 1/(1+r);              % Discount factor
gamma = 1.5;                 % Coefficient of relative risk aversion
startA = 1;                  % How much asset do people start life with

%% ------------------------------------------------------------------------ 
% GRIDS
% choose dimensiona and select methods to construct grids

%The grid for assets
%--------------------------------------%
numPointsA = 20;             % number of points in the discretised asset grid
gridMethod = 'equalsteps';   % method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps


%% ------------------------------------------------------------------------ 
% GET ASSET GRID
% populate grid for assets using 'gridMethod'
[ MinAss, MaxAss ] = getMinAndMaxAss(startA);
Agrid = NaN(T+1, numPointsA);
for ixt = 1:1:T+1
    Agrid(ixt, :) = getGrid(MinAss(ixt), MaxAss(ixt), numPointsA, gridMethod);
end

%% ------------------------------------------------------------------------ 
% SOLVE CONSUMER'S PROBLEM
% Get policy function and value function 

tic;        % start the clock
[ policyA1, policyC, val] = solveValueFunction;
toc

%% ------------------------------------------------------------------------
% Get the analytical policy functions if we have no uncertainty and
% borrowing is allowed

[ policyA1_analytic, policyC_analytic] = getPolicy_analytical;

% Get the ratio of numerical policy function to analytic policy functions
ratioOfPolicyC = policyC./policyC_analytic;
ratioOfPolicyA = policyA1./policyA1_analytic;

%% ------------------------------------------------------------------------ 
% SIMULATE CONSUMER'S PATHS
% start from initial level of assets and simulate optimal consumption and
% savings profiles over lifecycle

[ cpath, apath, vpath ] = simNoUncer(policyA1, val, startA);

%% ------------------------------------------------------------------------ 
% PLOTS
% Plots some features of the simulation and simulation

% Plot paths of consumption, income and assets 
plotNumber = 0;
plotCpath(cpath)
plotCZoomedOut(cpath)
plotApath(apath)

% Now plot value and policy functions
whichyear = 1;
plotNode1 = 3;
plotNodeLast = numPointsA; 
plots;

% ------------------------------------------------------------------------ 
% ------------------------------------------------------------------------ 
