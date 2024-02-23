## Life Cycle Model of Consumption and Savings with Two Assets

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
## Test Cases
################################################################################


# Test 1: do they only use the liquid asset?
# Note: must make transaction cost apply after retirement too
# model = setmodel(;beta = 0.95, r_b = 0.05, r = 0.05, adj_cost_fixed = 0.0, adj_cost_prop = 0.4, det_inc = false, rho = 0.75, adj_cost_post_ret = true)

# Test 2: do they only use the illiquid asset?
# model = setmodel(;beta = 0.95, r_b = 0.05, r = 0.0, adj_cost = false, det_inc = false, rho = 0.75)

# Test 3: they should put nothing in the illiquid account.... what's going on here !??
# model = setmodel(;beta = 0.95, r_b = 0.04, r = 0.05, adj_cost_fixed = 0.2, adj_cost_prop = 0.8, det_inc = false, rho = 0.2)

# Test 4: they should put nothing in the illiquid account when it's return dominated
# model = setmodel(;beta = 0.95, r_b = 0.02, r = 0.05, adj_cost_fixed = 0.0, adj_cost_prop = 0.0, det_inc = false, rho = 0.2)

# Test 1 w/ temptation: given equal returns, do they use ret account when tempted?
# model = setmodel(;beta = 0.95, lambda = 0.2, r_b = 0.05, r = 0.05, adj_cost_fixed = 0.0, adj_cost_prop = 0.4, det_inc = false, rho = 0.75, adj_cost_post_ret = true)
# model = setmodel(;beta = 0.95, lambda = 0.1, r_b = 0.05, r = 0.05, adj_cost_fixed = 0.0, adj_cost_prop = 0.6, det_inc = false)

################################################################################
# Full Model
################################################################################

# Now add hump shaped income profile
# model = setmodel(;beta = 0.95, r_b = 0.05, r = 0.0, adj_cost = false)

# Model w/out Temptation
model = setmodel(;beta = 0.95, gamma = 1.5, r_b = 0.04, r = 0.0)

# Model w/ Temptation
# model = setmodel(;beta = 0.95, lambda = 0.2, gamma = 1.5, r_b = 0.04, r = 0.0)

################################################################################
# Constants
################################################################################

## If test cases:
# const T                    = 10                  # Number of time period
# const Tretire              = 7                   # Age at which retirement happens

## If full model:
const T                    = 60                  # Number of time period
const Tretire              = 45                   # Age at which retirement happens
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

################################################################################
## Solve & Simulate
################################################################################

println("Solve Value Function: Parallel")
@time policyA1, policyB1, policyC, V, EV, V_NA  = solveValueFunctionPar(model, Agrid, Bgrid, Ygrid, incTransitionMrx)
# @time policyA1, policyB1, policyC, V, EV, V_NA  = solveValueFunctionPar(model, Agrid, Bgrid, Ygrid, incTransitionMrx)

@time cpath, apath, bpath, vpath, ypath, ewpath = simWithUncer(model, Agrid, Bgrid, Ygrid, policyA1, policyB1, EV)
# @time cpath, apath, bpath, vpath, ypath, ewpath = simWithUncer(model, Agrid, Bgrid, Ygrid, policyA1, policyB1, EV)

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


# Plot frequency of early withdrawals
ew_mean = mean(ewpath, dims =2)
plt = plot( [1:length(ewpath[:, 1])], ew_mean, linewidth = 2, label = "Early Withdrawal")
display(plt)


# Plot Life Cycle Profile for 3 HHs
for ixHH = 1:3 # [1, 994] #  994:994
    plt = plot([1:length(cpath[:,  1])], cpath[:, ixHH],linewidth = 2, label = "Consumption")
    plt = plot!([1:length(apath[:, 1])], apath[:, ixHH], linewidth = 2, label="Liq Assets")
    plt = plot!([1:length(bpath[:, 1])], bpath[:, ixHH], linewidth = 2, label="Illiq Assets")
    plt = plot!([1:length(ypath[:, 1])], ypath[:, ixHH],linewidth = 2, label = "Income", xlabel="Age", title="Household "*string(ixHH))

    # Shade early withdrawal
    # plt = plot!([1:length(ewpath[:, 1])], ewpath[:, ixHH],linewidth = 2, label = "Early Withdrawal")
    # plt = vspan!(collect(1:T)[ewpath[:, ixHH] .== 1.0]; alpha = 0.2, label = "Early Withdrawal")
    plt = vspan!( [[i-.5, i+.5] for i=1:T][ewpath[:, ixHH] .== 1.0]; alpha = 0.1, color = :gray, label ="") #, label = "Early Withdrawal")

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



