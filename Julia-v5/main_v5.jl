## Simple Life Cycle Model of Consumption and Savings
## Based off of matlab code by Monica Costa Dias and Cormac O'Dea

# Model v1: finite consumption saving endowment problem
# Model v3: add deterministic income stream
# Model v5: realistic income uncertainty

# This code was written in Julia v0.6.2 in September 2017.
# Updated for Julia v1.8 in April 2023.
# This code is designed to be easy to understand, not to be as fast as possible

################################################################################
## Load packages
################################################################################

using Distributed, Interpolations, QuadGK, Optim, Roots, LinearAlgebra, Random, Plots
import Random.seed!

################################################################################
## Run multi-threaded or not
################################################################################

# haven't actually run this with more than one thread 
runparallel = false

################################################################################
## Load dependencies
################################################################################

include("src/modelSetup.jl")
include("src/utils.jl")
include("src/model.jl")
include("src/modelEulerEquation.jl")
include("src/solveValueFunction.jl")
include("src/solveEulerEquation.jl")
include("src/solveValueFunctionPar.jl") # parallel version
include("src/simulation.jl")

################################################################################
## Define parameters and constants
################################################################################

# Define the parameters as a dictionary
# TODO: how does speed compare to using a constant?
# TODO: create a Dict of various objects (params, objs, etc)
# TODO: https://docs.julialang.org/en/stable/manual/performance-tips/#tools-1
params                     = Dict{String, Float64}()
params["tol"]              = 1e-5               # max allowed error
params["minCons"]          = 1e-5                # min allowed consumption
params["r"]                = 1.0/0.95 - 1.0      # Interest rate
params["beta"]             = 0.95                # 1/(1+r) # Discount factor
params["gamma"]            = 1.5                 # Coefficient of relative risk aversion
params["gamma_mod"]        = 1.0-params["gamma"] # For speed, just do this once
params["startA"]           = 0.0                 # How much asset do people start life with
params["mu"]               = 0.0                 # mean of initial log income
params["sigma"]            = 0.25                # variance of innovations to log income
params["rho"]              = 0.75                # persistency of log income

# Constants
const interpMethod         = "linear"            # for now, I only allow linear option
const T                    = 60                  # Number of time period
const Tretire              = 45                  # Age at which retirement happens
const borrowingAllowed     = 0                   # allow borrowing
const isUncertainty        = 1                   # uncertain income (currently: only works if isUncertainty == 1)
const numPointsY           = 9                   # number of points in the income grid
const numPointsA           = 100                  # number of points in the discretised asset grid
const gridMethod           = "5logsteps"         # method to construct grid. One of equalsteps or 5logsteps
const normBnd              = 3                   # truncate the normal distrib: ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma
const numSims              = 10                  # How many individuals to simulate
const useEulerEquation     = false               # Solve the model using the euler equation?
const saveValue_inEE       = true               # When using euler equation to solve the model, do we want to compute EV? (Note: adds time due to interpolation)
const linearise            = true               # Whether to linearise the slope of EdU when using EE


################################################################################
## Setup Model
################################################################################

# Get income grid
Ygrid, incTransitionMrx, minInc, maxInc = getIncomeGrid(params)

# Get asset grid
MinAss, MaxAss = getMinAndMaxAss(params, minInc, maxInc)
Agrid = zeros(T+1, numPointsA)
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
end

################################################################################
## Solve
################################################################################

if useEulerEquation == true
    println("Solve Euler Equation")
    #@time policyA1, policyC, V, EV, dU, EdU = solveEulerEquation(params, Agrid, Ygrid, incTransitionMrx)
    @time policyA1, policyC, V, EV, dU, EdU = solveEulerEquation(params, Agrid, Ygrid, incTransitionMrx) # NOTE: julia uses just in time compilation, so the second run is faster than the first because the code has now been compiled
elseif runparallel == false
    println("Solve Value Function: Serial")
    @time policyA1, policyC, V, EV  = solveValueFunction(params, Agrid, Ygrid, incTransitionMrx)
else
    println("Solve Value Function: Parallel")
    @time policyA1, policyC, V, EV  = solveValueFunctionPar(params, Agrid, Ygrid, incTransitionMrx)
end

# NOTE: evaluation time will be faster if you run it a second time, due to just in time compilation

################################################################################
## Simulate
################################################################################

cpath, apath, vpath, ypath = simWithUncer(params, Agrid, Ygrid, policyA1, EV)

################################################################################
## Plot
################################################################################

# Plot example
plot([1:length(cpath[:, 1])], cpath[:,1],linewidth = 2)
plot!([1:length(apath[:, 1])], apath[:,1], linewidth = 2)
plot!([1:length(ypath[:, 1])], ypath[:,1],linewidth = 2)

