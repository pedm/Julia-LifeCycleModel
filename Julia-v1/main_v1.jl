# v1: finite consumption saving endowment problem

using Interpolations
using Optim

include("src/modelSetup.jl")
include("src/model.jl")
include("src/solveValueFunction.jl")
include("src/simulation.jl")
include("src/plots.jl")

const minCons = 1e-5              # min allowed consumption
const T = 80                      # Number of time period
const r = 0.01                    # Interest rate
const beta = 1/(1+r)              # Discount factor
const gamma = 1.5                 # Coefficient of relative risk aversion
const gamma_mod = 1.0-gamma       # For speed, just do this once
const startA = 1                  # How much asset do people start life with
const numPointsA = 20             # number of points in the discretised asset grid
gridMethod = "equalsteps"         # method to construct grid. One of equalsteps or 5logsteps

# GET ASSET GRID
# populate grid for assets using 'gridMethod'
MinAss, MaxAss = getMinAndMaxAss()
Agrid = zeros(T+1, numPointsA)
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
end

@time policyA1, policyC, V  = solveValueFunction(minCons, Agrid)

cpath, apath, vpath = simNoUncer(T, r, Agrid, policyA1, V, startA);

# Load PyPlot module
# ENV["MPLBACKEND"]="qt4agg"
# using PyPlot

using Plots
plotly()

# Plots
plotCpath(cpath)
plotApath(apath, MinAss)
plotCAndApaths( cpath, apath );
