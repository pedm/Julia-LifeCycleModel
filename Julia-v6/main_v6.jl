## Simple Life Cycle Model of Consumption and Savings
## Based off of matlab code by Monica Costa-Dias and Cormac O'Dea

# v1: finite consumption saving endowment problem
# v3: add deterministic income stream
# v5: realistic income uncertainty with persistent income shocks
# v6: and transitory income shocks

################################################################################
# Euler Equation Update (Jan 2019)

# Next steps:
# 1) add the option to linearize the slope of EdU when using EE (as in Cormac's code) -- HD working on it
# 2) for a given individual, plot the euler equation over time. is it flat?
# 3) plot marginal utility of consumption by A1 (given some income)
# 4) for some arbitrary set of states, can we plot eulerforzero() across A1? then decompose this into the
# 5) modify income process so that you get less income when young. turn on borrowing, check that households consumption smooth by borrowing when young. turn off borrowing, check that euler equation is violated when young due to credit constraints

# Way down the road
# 1) create a higher interest rate for borrowing than lending, solve with EE method (will require two EEs, as explained in KV2014 appendix)
# 2) allow for a discrete choice in the model, with some sort of smoothing, and solve with euler equation (following Blundell, Costa-Dias, Meghir, and Shaw)

# note that consumption in retirement is much smoother when solving using euler equation
# why does consumption go up in retirement when 30 asset points?
# did tony smith suggest a faster way to find roots when using interpolation?

################################################################################
## Run parallel or not
################################################################################

runparallel = false

if runparallel == true
    # Add workers
    println("Add workers")
    if nprocs() == 1
        addprocs(7)
    end
    println(nprocs())
end

################################################################################
## Load dependencies
################################################################################

# Load packages
try
	println("Loading packages")
	@everywhere using Interpolations
	@everywhere using Optim
	using Roots
    using QuadGK
    using FastGaussQuadrature
catch
	println("Installing packages")
	Pkg.add("Interpolations")
	Pkg.add("Optim")
	Pkg.add("Roots")
    Pkg.add("QuadGK")
    Pkg.add("FastGaussQuadrature")

	@everywhere using Interpolations
	@everywhere using Optim
    using FastGaussQuadrature
	using Roots
    using QuadGK
    using FastGaussQuadrature
end

include("src/modelSetup.jl")
include("src/utils.jl")
include("src/model.jl")
include("src/modelEulerEquation.jl")
include("src/solveValueFunction.jl")
include("src/solveEulerEquation.jl")
include("src/solveValueFunctionPar.jl") # parallel version
include("src/simulation.jl")
include("src/plots.jl")

################################################################################
## Define parameters and constants
################################################################################

# Define the parameters as a dictionary
# TODO: how does speed compare to using a constant?
# TODO: create a Dict of various objects (params, objs, etc)
# TODO: https://docs.julialang.org/en/stable/manual/performance-tips/#tools-1
params                     = Dict{String, Float64}()
params["tol"]              = 1e-10               # max allowed error
params["minCons"]          = 1e-4                # min allowed consumption
params["r"]                = 1.0/0.95 - 1.0      # Interest rate
params["beta"]             = 0.95                # 1/(1+r) # Discount factor
params["gamma"]            = 2.0                 # Coefficient of relative risk aversion
params["gamma_mod"]        = 1.0-params["gamma"] # For speed, just do this once
params["startA"]           = 0.0                 # How much asset do people start life with
params["mu"]               = 0.0                 # mean of initial log income
params["sigma"]            = 0.25                # variance of innovations to log income
params["rho"]              = 0.95                # persistency of log income

# Seems something is not working right -- why does consumption go up so much with time???
params["sigma"]            = 0.3                # variance of innovations to log income
params["sigma_trans"]      = 0.000000001           # variance of transitory innovations to log income


params["sigma"]            = 0.000000001                # variance of innovations to log income
params["sigma_trans"]      = 0.000000001           # variance of transitory innovations to log income


# Constants
const interpMethod         = "linear"            # for now, I only allow linear option
const T                    = 60                  # Number of time period
const Tretire              = 45                  # Age at which retirement happens
const borrowingAllowed     = 0                   # allow borrowing
const isUncertainty        = 1                   # uncertain income (currently: only works if isUncertainty == 1)
const numPointsY           = 10                  # number of points in the income grid
const numPointsA           = 50                 # number of points in the discretised asset grid
const numPointsYTrans      = 15                   # number of points in the transitory income grid
const gridMethod           = "5logsteps"         # method to construct grid. One of equalsteps or 5logsteps
const normBnd              = 3                   # truncate the normal distrib: ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma
const numSims              = 10                  # How many individuals to simulate
const useEulerEquation     = true                # Solve the model using the euler equation?
const saveValue_inEE       = true                # When using euler equation to solve the model, do we want to compute EV? (Note: adds time due to interpolation)

################################################################################
## Setup Model
################################################################################

# Get income grid
Ygrid, incTransitionMrx, minInc, maxInc = getIncomeGrid(params)

μtransshocks = -0.5 * params["sigma_trans"] # this comes from Zeldes 1989 and Kovacs 2015. The expected value of a log normal variable with mean μ and variance σ2 is given by exp(μ + σ2/2)
Ytrans_grid, Ytrans_weights = gausshermite_normal_distribution(numPointsYTrans, μtransshocks, params["sigma_trans"] )
Ytrans_grid = exp.(Ytrans_grid)
maxInc = maxInc .* Ytrans_grid[end]

# Get asset grid
MinAss, MaxAss = getMinAndMaxAss(params, minInc, maxInc)
Agrid = zeros(T+1, numPointsA)
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
end
# QUESTION: why does the smallest value in Agrid not correspond to smallest value in MinAss?

################################################################################
## Solve
################################################################################

if useEulerEquation
    println("Solve Euler Equation")
    @time policyA1, policyC, V, EV, dU, EdU = solveEulerEquation(params, Agrid, Ygrid, incTransitionMrx)
    @time policyA1, policyC, V, EV, dU, EdU = solveEulerEquation(params, Agrid, Ygrid, incTransitionMrx) # NOTE: julia uses just in time compilation, so the second run is faster than the first because the code has now been compiled
elseif runparallel == false
    println("Solve Value Function: Serial")
    @time policyA1, policyC, V, EV  = solveValueFunction(params, Agrid, Ygrid, incTransitionMrx)
else
    println("Solve Value Function: Parallel")
    @time policyA1, policyC, V, EV  = solveValueFunctionPar(params, Agrid, Ygrid, incTransitionMrx)
end

## NOTE: evaluation time will be faster if you run it a second time, due to just in time compilation

################################################################################
## Simulate
################################################################################

cpath, apath, vpath, ypath = simWithUncer(params, Agrid, Ygrid, policyA1, EV)

################################################################################
## Plot
################################################################################

try
	using Plots
catch
	println("Installing packages")
	Pkg.add("Plots")
	using Plots
end

# http://docs.juliaplots.org/latest/
plotly() # this works better for development -- shows the plot in firefox. but it does not allow you to save as a pdf or png
# gr() # this works better on the server -- does not display the plot, but allows you to save it
# or try with plotlyjs()

plot_policyA1(10)
plot_policyA1(55, 1)
plot_V(5)
plot_EV_over_time()

plotCpath(cpath)
# savefig("TEST_cpath.pdf")

plotApath(apath, MinAss)
# savefig("TEST_apath.pdf")

# plotYAndCpaths( ypath, cpath );
# savefig("TEST_ycpath.pdf")

plotYCAndApaths( ypath, cpath, apath );

# BUG! When both sigmas set to 0.000000001, why does consumption go up in retirement???

################################################################################
## Profile
################################################################################

# Next steps
# TODO: get simNoUncer() working in this version
# TODO: dont use zeros(), instead use a = Array(Float64, size)
# TODO: use profiling:

## Tools to optimize code:
# @timev solveValueFunction(params, Agrid, Ygrid, incTransitionMrx);
# @code_warntype solveValueFunction(params, Agrid, Ygrid, incTransitionMrx);

# Profile.clear()  # in case we have any previous profiling data
# @profile solveValueFunction(params, Agrid, Ygrid, incTransitionMrx)
# Profile.print()

# try
# 	using ProfileView
# catch
# 	println("Installing packages")
# 	Pkg.add("ProfileView")
# 	using ProfileView
# end
#
# ProfileView.view()
