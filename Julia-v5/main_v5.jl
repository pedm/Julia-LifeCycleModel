# v1: finite consumption saving endowment problem
# v3: add deterministic income stream
# v5: realistic income uncertainty
# v5_improved: dictionary of params, better plots, get parallel working
#   March 6th Update - replace constants with a dictionary of params
#   This will be better down the road when we do estimation
#   March 21 Update - better plots

# TODO: work on a better parallel version

################################################################################
## Run parallel or not
################################################################################

runparallel = true

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
catch
	println("Installing packages")
	Pkg.add("Interpolations")
	Pkg.add("Optim")
	Pkg.add("Roots")
    Pkg.add("QuadGK")

	@everywhere using Interpolations
	@everywhere using Optim
	using Roots
    using QuadGK
end

include("src/modelSetup.jl")
include("src/model.jl")
include("src/solveValueFunction.jl")
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
params["minCons"]          = 1e-5                # min allowed consumption
params["r"]                = 0.02                # Interest rate
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
const numPointsY           = 100                  # number of points in the income grid
const numPointsA           = 100                  # number of points in the discretised asset grid
const gridMethod           = "5logsteps"         # method to construct grid. One of equalsteps or 5logsteps
const normBnd              = 3                   # truncate the normal distrib: ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma
const numSims              = 10                  # How many individuals to simulate

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

if runparallel == false
    println("First pass")
    @time policyA1, policyC, V, EV  = solveValueFunction(params, Agrid, Ygrid, incTransitionMrx)

    println("Second pass")
    @time policyA1, policyC, V, EV  = solveValueFunction(params, Agrid, Ygrid, incTransitionMrx)

else
    println("First pass")
    @time policyA1, policyC, V, EV  = solveValueFunctionPar(params, Agrid, Ygrid, incTransitionMrx)

    println("Second pass")
    @time policyA1, policyC, V, EV  = solveValueFunctionPar(params, Agrid, Ygrid, incTransitionMrx)
end

## Look at evaluation time

# @elapsed 1
# # Average time
# time = zeros(30)
# for i=1:30
#     time[i] = @elapsed policyA1, policyC, V, EV  = solveValueFunction(params, Agrid, Ygrid, incTransitionMrx);
# end
# println("Mean Execution Time")
# println(mean(time))

################################################################################
## Simulate
################################################################################

cpath, apath, vpath, ypath = simWithUncer(params, Agrid, Ygrid, policyA1, EV)

## Tools to optimize code:
# @timev solveValueFunction(params, Agrid, Ygrid, incTransitionMrx);
# @code_warntype solveValueFunction(params, Agrid, Ygrid, incTransitionMrx);

# TODO: why is negV ANY rather than Float64?

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

plotCpath(cpath)
# savefig("TEST_cpath.pdf")

plotApath(apath, MinAss)
# savefig("TEST_apath.pdf")

# plotYAndCpaths( ypath, cpath );
# savefig("TEST_ycpath.pdf")

plotYCAndApaths( ypath, cpath, apath );


# Next steps
# TODO: get simNoUncer() working in this version
# TODO: dont use zeros(), instead use a = Array(Float64, size)
# TODO: use profiling:

################################################################################
## Profile
################################################################################

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
