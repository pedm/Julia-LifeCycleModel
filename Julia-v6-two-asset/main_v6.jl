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

@time cpath, apath, bpath, vpath, ypath, ewpath = simWithUncer(model, Agrid, Bgrid, Ygrid, policyA1, policyB1, EV)

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

## Event Study:

# Find all cases where income falls by 20% or more 

pct_change_y = (ypath[2:Tretire-1, :] - ypath[1:Tretire-2, :]) ./ ypath[1:Tretire-2, :]
mean(pct_change_y .< -0.2)


# DataFrames Tutorial: https://github.com/bkamins/Julia-DataFrames-Tutorial/tree/master
# Tidier: 

using DataFrames 
using ShiftedArrays
using PanelDataTools

N = model["ints"]["numSims"]
id = repeat(collect(1:N)', T, 1)
t = repeat(collect(1:T), 1, N)

df = DataFrame(id = id[:], t = t[:], y = ypath[:], c = cpath[:], a = apath[1:T, :][:], b = bpath[1:T, :][:], early_w = ewpath[:] )

# Not needed / not working:
# groupby(df, :id)
# transform!(groupby(df, :id), :y => (x -> log.(x) .- log.(lag(x))) => :d_y)
# transform!(groupby(df, :id), :y => (x -> log.(x) .- log.(lag(x))) => :d_y)

# Use Eirik's package
paneldf!(df,:id,:t)
lag!(df,:y) # Create lag of income
df[!, :d_y] = df[!, :y]- df[!, :L1y]
df[!, :pct_change_y] = df[!, :d_y] ./ df[!, :L1y]

# Create event time: time to/from first fall in income of 20% or more
df[!, :event]      = df[!, :pct_change_y] .< -0.2
df[!, :event_flag] = 1(df[!, :event])
df[ismissing.(df.event_flag), :event_flag] .= 0


# # Select the first event_flag by id
# df2 = transform(groupby(df, :id)) do group
#     group.event_flag  = coalesce.(group.event_flag, 0) # Replace missing values with 0
#     group.first_event = ones(size(group.event_flag)) * findfirst(==(1), group.event_flag)
#     group.event_time  = group.t - Int64.(group.first_event)
#     return group
# end

# Ensure that they are not too close to retirement
df2 = filter(row -> 5 <= row.t <= Tretire - 10, df2)

# Modify the code to randomly select one occurrence of event_flag == 1 for each id
df2 = transform(groupby(df, :id)) do group
    # Find indices where event_flag == 1 
    # Also ensure at least 5 years away from retirement
    event_indices = findall((group.event_flag .== 1) )

    # Randomly select one index
    random_index = rand(1:length(event_indices))

    # Get the time of the randomly selected event
    event_time = group.t[event_indices[random_index]]
    # Calculate event_time relative to the randomly selected event
    group.event_time = group.t .- event_time
    return group
end



filtered_df = filter(row -> -5 <= row.event_time <= 5, df2)

# Make a balanced panel: only keep the ids that have all observations for time -5 to 5
grouped_df = groupby(filtered_df, :id)                                            # Group the filtered DataFrame by the id column
filtered_groups = filter(grp -> length(unique(grp.event_time)) == 11, grouped_df) # Filter groups to keep only those with observations for time periods -5 to 5
final_df = vcat(filtered_groups...)                                               # Combine the filtered groups back into a DataFrame


# Collapse by time and compute the mean of variable x
collapsed_df = combine(groupby(final_df, :event_time)) do group
    DataFrame(y = mean(group.y), pct_change_y = mean(group.pct_change_y), age = mean(group.t) , c = mean(group.c), a = mean(group.a), b = mean(group.b), early_w = mean(group.early_w))
end


# Create plots:
# plot(collapsed_df.event_time, collapsed_df.age)

plot(collapsed_df.event_time, collapsed_df.age)
plot(collapsed_df.event_time, collapsed_df.y)

p1 = plot(collapsed_df.event_time, collapsed_df.pct_change_y, title = "Change Income")
p2 = plot(collapsed_df.event_time, collapsed_df.c, title = "Consumption")

p3 = plot(collapsed_df.event_time, collapsed_df.a, title = "Liquid Assets")
p4 = plot(collapsed_df.event_time, collapsed_df.b, title = "Retirement Assets")

p5 = plot(collapsed_df.event_time, collapsed_df.y .- collapsed_df.c, title = "Saving Rate (Y-C)", xlabel = "Event Time (Neg Shock at 0)")
p6 = plot(collapsed_df.event_time, collapsed_df.early_w, title = "Early Withdrawal Probability", xlabel = "Event Time (Neg Shock at 0)")

# p4 = plot(collapsed_df.event_time, collapsed_df.a ./ collapsed_df.y, title = "Fin Assets / Y")


l = @layout [a b; c d; e f]
plot(p1, p2, p3, p4, p5, p6, layout = l, size = (900, 800), legend = false)


#############################################################################################
# Now create event study with heterogeneity by liq assets / income ratio at time t=-1

filtered_df = filter(row -> row.event_time == -1, final_df) # Find all obs at t = -1
filtered_df.a_over_y = filtered_df.a ./ filtered_df.y # Compute assets over income

# Step 1a: Assign to quartiles - unconditional on age
# q = quantile(filtered_df.a_over_y, [0.25, 0.5, 0.75])
# filtered_df[!, :quartile_dummy] .= ifelse.(filtered_df.event_time .== -1,
#                                         ifelse.(filtered_df.a_over_y .<= q[1], 1,
#                                         ifelse.(filtered_df.a_over_y .<= q[2], 2,
#                                         ifelse.(filtered_df.a_over_y .<= q[3], 3, 4))),
#                                         missing)


# Step 1b: Assign to quartiles - conditional on age
quartiles_by_age = combine(groupby(filtered_df, :t)) do group
    q = quantile(group.a_over_y, [0.25, 0.5, 0.75])
    DataFrame(t = group.t[1], Q1 = q[1], Q2 = q[2], Q3 = q[3])
end
filtered_df = leftjoin(filtered_df, quartiles_by_age, on = :t)
filtered_df[!, :quartile_dummy] .= ifelse.(filtered_df.event_time .== -1,
                                        ifelse.(filtered_df.a_over_y .<= filtered_df.Q1, 1,
                                        ifelse.(filtered_df.a_over_y .<= filtered_df.Q2, 2,
                                        ifelse.(filtered_df.a_over_y .<= filtered_df.Q3, 3, 4))),
                                        missing)

# Step 2: Now just select the columns I want
quartiles_by_id = select(filtered_df, [:id, :quartile_dummy])

# Step 3: Merge quartiles back to the original DataFrame based on the id
merged_df = leftjoin(final_df, quartiles_by_id, on = :id)



# Collapse by time and compute the mean of variable x
collapsed_df = combine(groupby(merged_df, [:event_time, :quartile_dummy] )) do group
    DataFrame(y = mean(group.y), pct_change_y = mean(group.pct_change_y), age = mean(group.t) , c = mean(group.c), a = mean(group.a), b = mean(group.b), early_w = mean(group.early_w))
end


# Create plots:
# plot(collapsed_df.event_time, collapsed_df.age)

# Get unique IDs and time periods
times = unique(collapsed_df.event_time)
quart = unique(collapsed_df.quartile_dummy)

# Reshape DataFrame to an array
time_array = [collapsed_df[(collapsed_df.quartile_dummy .== q) .& (collapsed_df.event_time .== t), :].event_time[1] for t in times, q in quart]

get_array(var) = [collapsed_df[(collapsed_df.quartile_dummy .== q) .& (collapsed_df.event_time .== t), :][!, var][1] for t in times, q in quart]


p1 = plot(time_array, get_array(:pct_change_y), title = "Change Income")
p2 = plot(time_array, get_array(:c), title = "Consumption")
p3 = plot(time_array, get_array(:a), title = "Liquid Assets")
p4 = plot(time_array, get_array(:b), title = "Retirement Assets")
p5 = plot(time_array, get_array(:y) .- get_array(:c), title = "Saving Rate (Y-C)", xlabel = "Event Time (Neg Shock at 0)")
p6 = plot(time_array, get_array(:early_w), title = "Early Withdrawal Probability", xlabel = "Event Time (Neg Shock at 0)")

p7 = plot(time_array, get_array(:age), title = "Age", label = ["Q1" "Q2" "Q3" "Q4"])
p8 = plot(time_array, get_array(:y),   title = "Income", label = ["Q1" "Q2" "Q3" "Q4"])


l = @layout [a b; c d; e f; g h]
plot(p1, p2, p3, p4, p7, p8, p5, p6, layout = l, size = (800, 900), label = ["Q1" "Q2" "Q3" "Q4"])

l = @layout [a b; c d; e f]
plot(p1, p2, p3, p4, p5, p6, layout = l, size = (800, 900), label = ["Q1" "Q2" "Q3" "Q4"])


################################################################################################################

## Normalize to event time t = -1 
# Filter to find the value of variable x at time t = -1 within each quartile
x_at_t_minus_1 = combine(groupby(collapsed_df[collapsed_df.event_time .== -1, :], :quartile_dummy)) do group
    DataFrame(quartile_dummy = group.quartile_dummy[1], c_at_t_minus_1 = group.c)
end

# Merge x_at_t_minus_1 back to the original DataFrame based on the quartile
merged_df = leftjoin(collapsed_df, x_at_t_minus_1, on = :quartile_dummy)

# Normalize variable x based on its value at time t = -1 within each quartile
merged_df[!, :normalized_c] = (merged_df.c .- merged_df.c_at_t_minus_1) ./ merged_df.c_at_t_minus_1

# Display the updated DataFrame
show(merged_df)


get_array_normalized(var) = [merged_df[(merged_df.quartile_dummy .== q) .& (merged_df.event_time .== t), :][!, var][1] for t in times, q in quart]

p2 = plot(time_array, get_array_normalized(:normalized_c), title = "Normalized Consumption", label = ["Q1" "Q2" "Q3" "Q4"])
