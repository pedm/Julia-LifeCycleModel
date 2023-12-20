## Simple Life Cycle Model of Consumption and Savings

# Model v1: finite consumption saving endowment problem
# Model v3: add deterministic income stream
# Model v5: realistic income uncertainty
# Model v6: liquid and illiquid assets

# This code was written in Julia v0.6.2 in September 2017.
# Updated for Julia v1.8 in April 2023.
# This code is designed to be easy to understand, not to be as fast as possible

################################################################################
## Load packages
################################################################################

using Distributed, Interpolations, QuadGK, Optim, Roots, LinearAlgebra, Random, Plots, Statistics
import Random.seed!

################################################################################
## Run multi-threaded or not
################################################################################

# true or false
runparallel = true

################################################################################
## Load dependencies
################################################################################

include("src/convenience.jl") 
include("src/modelSetup.jl")
include("src/utils.jl")
include("src/model.jl")
include("src/modelEulerEquation.jl")
include("src/solveValueFunctionPar.jl") # parallel version
include("src/simulation.jl")
include("src/rouwenhorst/rouwenhorst.jl")


################################################################################
## Define parameters and constants
################################################################################

# Test 1: do they only use the liquid asset?
# Note: also helps to make transaction cost apply after retirement too
model = setmodel(;beta = 0.95, r_b = 0.05, r = 0.05, adj_cost_fixed = 0.0, adj_cost_prop = 0.4, det_inc = false, rho = 0.75)

# Test 2: do they only use the illiquid asset?
# model = setmodel(;beta = 0.95, r_b = 0.05, r = 0.0, adj_cost = false, det_inc = false, rho = 0.75)

# Now add hump shaped income profile
# model = setmodel(;beta = 0.95, r_b = 0.05, r = 0.0, adj_cost = false)

# Baseline: do they use some mix of both?
# model = setmodel(;beta = 0.95, gamma = 2, r_b = 0.04, r = 0.0)

# Constants
const T                    = 10                  # Number of time period
const Tretire              = 7                   # Age at which retirement happens
const numPointsY           = 5                   # number of points in the income grid
const numPointsA           = 60                  # number of points in the discretised asset grid -- seems helpful to have more liquid points, since it's used in the intermediate step
const numPointsB           = 50                  # number of points in the discretised asset grid


################################################################################
## Setup Model
################################################################################

# Get income grid
Ygrid, incTransitionMrx, minInc, maxInc, det_income = getIncomeGrid(model)

# Get asset grids
Agrid, Bgrid = getAssetGrid(model)

# for k in 1:3, j in 1:3, i in 1:3
#     @show (i, j, k)
# end

################################################################################
## Solve
################################################################################

println("Solve Value Function: Parallel")
@time policyA1, policyB1, policyC, V, EV, V_NA  = solveValueFunctionPar(model, Agrid, Bgrid, Ygrid, incTransitionMrx)
# @time policyA1, policyB1, policyC, V, EV, V_NA  = solveValueFunctionPar(model, Agrid, Bgrid, Ygrid, incTransitionMrx)


################################################################################
## Simulate
################################################################################

@time cpath, apath, bpath, vpath, ypath = simWithUncer(model, Agrid, Bgrid, Ygrid, policyA1, policyB1, EV)
@time cpath, apath, bpath, vpath, ypath = simWithUncer(model, Agrid, Bgrid, Ygrid, policyA1, policyB1, EV)



################################################################################
# LATER: Tests
################################################################################

# For Test #1: this should be zero everywhere
[maximum(policyB1[t, :, 1, :]) for t = 1:T ]
# [maximum(policyB1[t, 1:40, 1, :]) for t = 1:T] # at least it's zero when ppl have very little

# Note that the above works when there's no fixed cost... but fails when there is a fixed cost (due to non-convexities)

# Alt Test 1: do they only use the liquid asset? (more difficult b/c the fixed cost introduces non-convexities, which )
# model = setpar(;beta = 0.95, r_b = 0.04, r = 0.04, adj_cost_fixed = 0.1)

# For Test 2:
# [maximum(policyA1[t, 1, 1:numPointsB-10, :]) for t = 1:T]

################################################################################
## Plot
################################################################################


# Plot Policy Fcn
ixt = T-1
# ixt = 2
ixY = 3
plot(Agrid[ixt,:], policyA1[ixt, :, 1:5:numPointsB, ixY], ylabel="Policy A1", xlabel="A0") # weird: illiq assets has no affect on A1 - so clearly the issue is in solution not sim
# plot(Agrid[ixt,:], policyA1[ixt, :, 5, :], ylabel="Policy A1", xlabel="A0") # in contrast, income has a clear effect on A1


### POLICY FCNS FOR B1 based on A0
plot(Agrid[ixt,:], policyB1[ixt, :, 1, ixY],  ylabel="Policy B1", xlabel="A0", title="Policy Fcn B1")


### POLICY FCNS FOR B1 based on B0
plot(Bgrid[ixt,:], policyB1[ixt, 1, :, ixY],  ylabel="Policy B1", xlabel="B0", title="Policy Fcn B1",     label="ixA0 = 1")
plot!(Bgrid[ixt,:], policyB1[ixt, 10, :, ixY], ylabel="Policy B1", xlabel="B0", title="Policy Fcn B1",    label="ixA0 = 10") 
plot!(Bgrid[ixt,:], policyB1[ixt, end-1, :, ixY], ylabel="Policy B1", xlabel="B0", title="Policy Fcn B1", label="ixA0 = end-1") 
plot!(Bgrid[ixt,:], policyB1[ixt, end, :, ixY], ylabel="Policy B1", xlabel="B0", title="Policy Fcn B1",   label="ixA0 = end") 
# Weird that for the household with max liquid assets, they set B1 = 0... WHY?



# plot(Agrid[ixt,:], policyAdj[ixt, :, 1:5:numPointsB, ixY], ylabel="Policy Adjust", xlabel="A0")



### EV PLOTS 
# Looks good: more liquid assets is a good thing
plot(Agrid[ixt,:], EV[ixt, :, 1:5:numPointsB, ixY], ylabel="EV", xlabel="A0")
plot(Agrid[ixt,:], EV[ixt, :, 1, ixY],              ylabel="EV", xlabel="A0")
plot(Bgrid[ixt,:], EV[ixt, 1, :, ixY],             ylabel="EV", xlabel="B0")

# Compare EV of A and B
ixt = 5
plot(Agrid[ixt,:], EV[ixt, :, 1, ixY],              ylabel="EV", label = "EV given A0 (and B0 = 0)", xlabel="A0")
plot!(Bgrid[ixt,:], EV[ixt, 1, :, ixY],             ylabel="EV", label = "EV given B0 (and A0 = 0)", xlabel="A0 / B0")


### V_NA PLOTS conditinal on Astar
ixt = 3
plot(1:length(V_NA[ixt, :, 1, ixY]), V_NA[ixt, :, 1, ixY],      ylabel="V_NA", xlabel="Astar", label="ixB0 = 1")
plot!(1:length(V_NA[ixt, :, 5, ixY]), V_NA[ixt, :, 1, ixY],     ylabel="V_NA", xlabel="Astar", label="ixB0 = 5")
plot!(1:length(V_NA[ixt, :, 10, ixY]), V_NA[ixt, :, 1, ixY],    ylabel="V_NA", xlabel="Astar", label="ixB0 = 10")
plot!(1:length(V_NA[ixt, :, end, ixY]), V_NA[ixt, :, 1, ixY],   ylabel="V_NA", xlabel="Astar", label="ixB0 = end")


### V_NA PLOTS conditinal on Bgrid - aka illiquid assets assuming you don't adjust
ixt = 3
#plt = plot(Bgrid[ixt,:], V_NA[ixt, 1, :, ixY], ylabel="V_NA", xlabel="B0", label = "ixAstar = 1")
plt = plot(Bgrid[ixt,:], V_NA[ixt, 2, :, ixY], ylabel="V_NA", xlabel="B0", label = "ixAstar = 2")
plt = plot!(Bgrid[ixt,:], V_NA[ixt, 5, :, ixY], ylabel="V_NA", xlabel="B0", label = "ixAstar = 5")
plt = plot!(Bgrid[ixt,:], V_NA[ixt, 10, :, ixY], ylabel="V_NA", xlabel="B0", label = "ixAstar = 10")
display(plt)
# seems this creates a weird effect where by far the highest V_NA is if you have very low B0.... why??
# Normal looking if ixAstar = 1, aka negative Astar
# TODO: why?
# Seems to happen if max_contrib = 40.... but not if max_contrib=1. Strange! 

# could keep trying to remove the extrapolation option from the solution method? dunno...

# SHIT! Still opposite direction of what we would expect....
# TODO: why does that go in the opposite direction?

# surf(collect(1:length(V_NA[ixt, :, 1, ixY])), Bgrid[ixt,:], V_NA[ixt, :, :, ixY],              ylabel="V_NA", xlabel="A0")

plot( maximum(bpath, dims = 2))
plot!(minimum(bpath, dims = 2))


plot( maximum(apath, dims = 2))
plot!(minimum(apath, dims = 2))

# TODO: look at value function given B... does it have curviture?

# Note that hhs start off life with B = 1... so no wonder it stays fixed at that level forever.
# Seems issue in solution method, since policy fcns make no sense

# Plot Life Cycle Profile for 3 HHs
for ixHH = 1:3
    plt = plot([1:length(cpath[:,  1])], cpath[:, ixHH],linewidth = 2, label = "Consumption")
    plt = plot!([1:length(apath[:, 1])], apath[:, ixHH], linewidth = 2, label="Liq Assets")
    plt = plot!([1:length(bpath[:, 1])], bpath[:, ixHH], linewidth = 2, label="Illiq Assets")
    plt = plot!([1:length(ypath[:, 1])], ypath[:, ixHH],linewidth = 2, label = "Income", xlabel="Age", title="Household "*string(ixHH))
    display(plt)
end

# Plot Average Life Cycle Profiles
cpath_mean = mean(cpath, dims = 2)
apath_mean = mean(apath, dims = 2)
bpath_mean = mean(bpath, dims = 2)
ypath_mean = mean(ypath, dims = 2)
plot( [1:length(cpath[:, 1])], cpath_mean, linewidth = 2, label = "Consumption")
plot!([1:length(apath[:, 1])], apath_mean, linewidth = 2, label = "Liq Assets")
plot!([1:length(bpath[:, 1])], bpath_mean, linewidth = 2, label = "Illiq Assets")
plot!([1:length(ypath[:, 1])], ypath_mean, linewidth = 2, label = "Income", xlabel="Age", title = "Life-Cycle Profile (Average)")



# Things to test:
# 1. If r = r_b, and there are adj costs, then ppl should only save in liquid
# 2. If no adjustment costs, and r_b > r, then ppl should only save in illiquid


# Q: is the issue the sim or solution? 
# Whats up with this max in bpath for so many households? it seems to level out at some point... odd
# how often does it extrapolate?
# Q: why is minimum(bpath) negative !? 

# DONE - Maybe I could try bounding y within the Ygrid - then turn off extrapolation in simulation method - seems issue not caused by this!
# Maybe I should turn off extrapolation in the solution?

# TODO: why does value of holding B conditional on not adjusting go in the opposite direction of what i expect? or does this make sense? dunno 
# TODO: perhaps put a cap on B1 holdings -- ensure that we're always on the A1star grid ?

# NEXT THINGS TO TRY:
# Maybe my issue comes from not keeping track of returns properly... seems sometimes i use end of period returns, othertimes start of period returns. Weird. 
# Also wait... why would B0 affect your spending at all in the "no adjust" case? Wouldn't it be better if it just drops out entirely? Yes I think so.... gotta try that. 

# Note: I have turned off "free withdrawals" after retirement -- now nobody should ever invest in the illiquid account. Weird!!

# What if i give income in the outer loop? 

# todo: 
# use policyA1_NA() and impose that B1 = 0.0 always... are they able to smooth consumption?

# TODO: in the simulation, could use the policy fcn for A1 conditional on B1 and atilde(B1)... would that be better? maybe?

# NOTES:
# Might want to play around with the gridMethod for Atildegrid. I am guessing we do not want so many points at zero, since it will be very rare for consumption to be zero. 



