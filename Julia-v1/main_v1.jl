## Simple Life Cycle Model of Consumption and Savings
## Based off of matlab code by Monica Costa Dias and Cormac O'Dea

# Model v1: finite consumption saving endowment problem
# Model v3: add deterministic income stream
# Model v5: realistic income uncertainty

# This code was written in Julia v0.6.2 in September 2017.
# Updates in Julia v1.8 in April 2023.
# This code is designed to be easy to understand, not to be as fast as possible

using Interpolations, Optim, Plots

include("src/modelSetup.jl")
include("src/model.jl")
include("src/solveValueFunction.jl")
include("src/simulation.jl")

# Define parameters
const minCons = 1e-5              # min allowed consumption
const T = 80                      # Number of time period
const r = 0.05                    # Interest rate
const beta = 1/(1+r)              # Discount factor
const gamma = 1.5                 # Coefficient of relative risk aversion
const gamma_mod = 1.0-gamma       # For speed, just do this once
const startA = 1                  # How much asset do people start life with
const numPointsA = 50            # number of points in the discretised asset grid
gridMethod = "5logsteps"         # method to construct grid. One of equalsteps, logsteps or 5logsteps

# GET ASSET GRID
# populate grid for assets using 'gridMethod'
MinAss, MaxAss = getMinAndMaxAss()
Agrid = zeros(T+1, numPointsA)
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
end

# Solve
policyA1, policyC, V  = solveValueFunction(minCons, Agrid)

# Simulate
cpath, apath, vpath = simNoUncer(T, r, Agrid, policyA1, V, startA);

# Plot
plot([1:length(cpath[:, 1])], cpath[:,1],linewidth = 2)
plot!([1:length(apath[:, 1])], apath[:,1], linewidth = 2)