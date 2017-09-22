# v1: finite consumption saving endowment problem
# v3: add deterministic income stream
# v5: realistic income uncertainty

using Interpolations
using Optim
using Roots

include("src/modelSetup.jl")
include("src/model.jl")
include("src/solveValueFunction.jl")
include("src/simulation.jl")
include("src/plots.jl")

const interpMethod = "linear"     # for now, I only allow linear option
const tol = 1e-10                 # max allowed error
const minCons = 1e-5              # min allowed consumption
const T = 60                      # Number of time period
const Tretire = 45                # Age at which retirement happens
const r = 0.02                    # Interest rate
const beta = 0.95 # 1/(1+r)       # Discount factor
const gamma = 1.5                 # Coefficient of relative risk aversion
const gamma_mod = 1.0-gamma       # For speed, just do this once
const startA = 0.0                # How much asset do people start life with
const mu = 0.0                    # mean of initial log income
const sigma = 0.25                # variance of innovations to log income
const rho = 0.75                  # persistency of log income
const borrowingAllowed = 0        # allow borrowing
const isUncertainty = 1           # uncertain income (currently: only works if isUncertainty == 1)
const numPointsY = 5              # number of points in the income grid
const numPointsA = 50             # number of points in the discretised asset grid
gridMethod = "5logsteps"          # method to construct grid. One of equalsteps or 5logsteps
const normBnd = 3                 # truncate the normal distrib: ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma
const numSims = 10                # How many individuals to simulate

# Get income grid
Ygrid, incTransitionMrx, minInc, maxInc = getIncomeGrid()

# Get asset grid
MinAss, MaxAss = getMinAndMaxAss(minInc, maxInc)
Agrid = zeros(T+1, numPointsA)
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
end

@time policyA1, policyC, V, EV  = solveValueFunction(Agrid, Ygrid, incTransitionMrx)

cpath, apath, vpath, ypath = simWithUncer(Agrid, Ygrid, policyA1, EV, startA)

ENV["MPLBACKEND"]="qt4agg"
using PyPlot
# PyPlot.close()

plotCpath(cpath)
plotApath(apath, MinAss)
# plotYAndCpaths( ypath, cpath );
plotYCAndApaths( ypath, cpath, apath );

# Next steps
# TODO: get simNoUncer() working in this version
# TODO: dont use zeros(), instead use a = Array(Float64, size)
# TODO: use profiling:

# Profile.clear()  # in case we have any previous profiling data
# @profile solveValueFunction(...)
# using ProfileView
# ProfileView.view()
