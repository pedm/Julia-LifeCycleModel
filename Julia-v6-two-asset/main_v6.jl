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

# true or fale
runparallel = true

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
params["r_b"]              = 1.0/0.98 - 1.0      # Interest rate
# params["r_b"]              = 0.05
params["r"]                = 0.01      # Interest rate

params["beta"]             = 0.98                # 1/(1+r) # Discount factor
params["gamma"]            = 1.5                 # Coefficient of relative risk aversion
params["gamma_mod"]        = 1.0-params["gamma"] # For speed, just do this once
params["startA"]           = 0.0                 # How much asset do people start life with
params["mu"]               = 0.0                 # mean of initial log income
params["sigma"]            = 0.25                # variance of innovations to log income
# params["sigma"]            = 0.05                # variance of innovations to log income
params["rho"]              = 0.75                # persistency of log income
params["adj_cost_fixed"]   = 0.1                 # fixed cost to adjust the illiquid asset - about 5% of avg annual income as the fixed cost
# params["adj_cost_prop"]    = 0.1                 # proportional cost to adjust the illiquid asset
params["adj_cost_prop"]    = 0.5                 # proportional cost to adjust the illiquid asset
params["max_contrib"]      = 0.5                 # maximum contribution to retirement account each period (arbitrary)
params["Yretire"]          = 0.1

# Constants
const interpMethod         = "linear"            # for now, I only allow linear option
const T                    = 60                  # Number of time period
const Tretire              = 45                  # Age at which retirement happens
# const T                    = 6                  # Number of time period
# const Tretire              = 4                  # Age at which retirement happensconst borrowingAllowed     = 0                   # allow borrowing
const isUncertainty        = 1                   # uncertain income (currently: only works if isUncertainty == 1)
const numPointsY           = 3                   # number of points in the income grid
const numPointsA           = 50                  # number of points in the discretised asset grid
const numPointsB           = 40                  # number of points in the discretised asset grid
const gridMethod           = "5logsteps"         # method to construct grid. One of equalsteps or 5logsteps
const normBnd              = 3                   # truncate the normal distrib: ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma
const numSims              = 500                  # How many individuals to simulate
const useEulerEquation     = false               # Solve the model using the euler equation?
const saveValue_inEE       = true               # When using euler equation to solve the model, do we want to compute EV? (Note: adds time due to interpolation)
const linearise            = true               # Whether to linearise the slope of EdU when using EE


################################################################################
## Setup Model
################################################################################

# Get income grid
Ygrid, incTransitionMrx, minInc, maxInc = getIncomeGrid(params)

# Get asset grids
MinAss, MaxAss, MaxB = getMinAndMaxAss(params, minInc, maxInc)
Agrid = zeros(T+1, numPointsA)
Bgrid = zeros(T+1, numPointsB)
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
    Bgrid[ixt, :] = getGrid(MinAss[ixt], MaxB[ixt], numPointsB, gridMethod)
end



# for k in 1:3, j in 1:3, i in 1:3
#     @show (i, j, k)
# end

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
    @time policyA1, policyB1, policyC, V, EV, V_NA, policyAdj  = solveValueFunctionPar(params, Agrid, Bgrid, Ygrid, incTransitionMrx)
end

# NOTE: evaluation time will be faster if you run it a second time, due to just in time compilation

################################################################################
## Simulate
################################################################################

cpath, apath, bpath, vpath, ypath = simWithUncer(params, Agrid, Bgrid, Ygrid, policyA1, policyB1, EV)

################################################################################
## Plot
################################################################################


# Plot Policy Fcn
ixt = T-1
# ixt = 2
# ixY = 3
plot(Agrid[ixt,:], policyA1[ixt, :, 1:5:numPointsB, ixY], ylabel="Policy A1", xlabel="A0") # weird: illiq assets has no affect on A1 - so clearly the issue is in solution not sim
# plot(Agrid[ixt,:], policyA1[ixt, :, 5, :], ylabel="Policy A1", xlabel="A0") # in contrast, income has a clear effect on A1


### POLICY FCNS FOR B1
plot(Agrid[ixt,:], policyB1[ixt, :, 1, ixY],  ylabel="Policy B1", xlabel="A0", title="Policy Fcn B1")
plot(Bgrid[ixt,:], policyB1[ixt, 1, :, ixY],  ylabel="Policy B1", xlabel="B0", title="Policy Fcn B1")
plot(Bgrid[ixt,:], policyB1[ixt, 10, :, ixY], ylabel="Policy B1", xlabel="B0", title="Policy Fcn B1") 

# plot(Agrid[ixt,:], policyAdj[ixt, :, 1:5:numPointsB, ixY], ylabel="Policy Adjust", xlabel="A0")



### EV PLOTS 
# Looks good: more liquid assets is a good thing
plot(Agrid[ixt,:], EV[ixt, :, 1:5:numPointsB, ixY], ylabel="EV", xlabel="A0")
plot(Agrid[ixt,:], EV[ixt, :, 1, ixY],              ylabel="EV", xlabel="A0")
plot(Bgrid[ixt,:], EV[ixt, 1, :, ixY],             ylabel="EV", xlabel="B0")

### V_NA PLOTS
plot(1:length(V_NA[ixt, :, 1, ixY]), V_NA[ixt, :, 1, ixY],              ylabel="V_NA", xlabel="A0")
plot(Bgrid[ixt,:], V_NA[ixt, 1, :, ixY],             ylabel="V_NA", xlabel="B0")

# SHIT! Still opposite direction of what we would expect....


# surf(collect(1:length(V_NA[ixt, :, 1, ixY])), Bgrid[ixt,:], V_NA[ixt, :, :, ixY],              ylabel="V_NA", xlabel="A0")


# TODO: look at value function given B... does it have curviture?

# Note that hhs start off life with B = 1... so no wonder it stays fixed at that level forever.
# Seems issue in solution method, since policy fcns make no sense

# Plot Life Cycle Profile for Indiv 1
plot([1:length(cpath[:, 1])], cpath[:,1],linewidth = 2, label = "Consumption")
plot!([1:length(apath[:, 1])], apath[:,1], linewidth = 2, label="Liq Assets")
plot!([1:length(bpath[:, 1])], bpath[:,1], linewidth = 2, label="Illiq Assets")
plot!([1:length(ypath[:, 1])], ypath[:,1],linewidth = 2, label = "Income", xlabel="Age")


# Plot Life Cycle Profile for Indiv 1
plot([1:length(cpath[:,  1])], cpath[:,2],linewidth = 2, label = "Consumption")
plot!([1:length(apath[:, 1])], apath[:,2], linewidth = 2, label="Liq Assets")
plot!([1:length(bpath[:, 1])], bpath[:,2], linewidth = 2, label="Illiq Assets")
plot!([1:length(ypath[:, 1])], ypath[:,2],linewidth = 2, label = "Income", xlabel="Age", title = "Life-Cycle Profile (HH 1)")


cpath_mean = mean(cpath, dims = 2)
apath_mean = mean(apath, dims = 2)
bpath_mean = mean(bpath, dims = 2)
ypath_mean = mean(ypath, dims = 2)
plot([1:length(cpath[:, 1])],  cpath_mean,linewidth = 2, label = "Consumption")
plot!([1:length(apath[:, 1])], apath_mean, linewidth = 2, label="Liq Assets")
plot!([1:length(bpath[:, 1])], bpath_mean, linewidth = 2, label="Illiq Assets")
plot!([1:length(ypath[:, 1])], ypath_mean,linewidth = 2, label = "Income", xlabel="Age", title = "Life-Cycle Profile (Avg)")
