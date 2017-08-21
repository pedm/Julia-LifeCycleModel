
using Interpolations
using Optim

interpMethod = "linear"     # for now, I only allow linear option
tol = 1e-10                 # max allowed error
minCons = 1e-5              # min allowed consumption
T = 80                      # Number of time period
r = 0.01                    # Interest rate
beta = 1/(1+r)              # Discount factor
startA = 1                  # How much asset do people start life with
numPointsA = 20             # number of points in the discretised asset grid
gridMethod = "equalsteps"   # method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps


function getMinAndMaxAss(T, r, minCons, startA)
    ## ------------------------------------------------------------------------
    # Initialise the output matrices
       BC   = zeros(T+1,1)
       maxA = zeros(T+1,1)
    # do I have to initialise the size NaN(T+1, 1) ?

    ## ------------------------------------------------------------------------
    # Iteratively, calculate the borrowing constraints, and maximum on asset
    # grid

    # Borrowing constraints
    BC[T + 1] = 0
    for ixt = T:-1:1
        BC[ixt] = BC[ixt+1]/(1+r) + minCons
    end

    # Maximum Assets
    maxA[1] = startA
    for ixt = 2:1:T+1
        maxA[ixt] = (maxA[ixt - 1] ) * (1+r)
    end
    return BC, maxA
end

function getGrid(minongrid, maxongrid, GridPoints, method)
    span = maxongrid - minongrid
    grid= linspace(minongrid, span, GridPoints)
end

# GET ASSET GRID
# populate grid for assets using 'gridMethod'
MinAss, MaxAss = getMinAndMaxAss(T, r, minCons, startA)
Agrid = zeros(T+1, numPointsA)
for ixt = 1:1:T+1
    Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
end


function utility(cons)
    gamma = 1.5                 # Coefficient of relative risk aversion
    # if gamma == 1
    #     utils = log(cons)
    # else
        utils = ((cons)^(1-gamma)  )/(1-gamma)
    # end
end


function objectivefunc(beta, r, itp, A1, A0)

    #-------------------------------------------------------------------------------#
    # This function returns the following quantity:
    # - (u(c) +  b V( A1))
    # where c is calculated from today's assets and tomorrow's assets

    #Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    #total value (u(c) + b * VA1

    cons = A0  - (A1)/(1+r)
    value = utility(cons) + beta * itp[A1]

    ## ------------------------------------------------------------------------
    #The optimisation routine that we will use searches for the minimum of the
    #function. We want the maximum. So we multiply out function here by -1 so
    #that the optimiser will fill the minimum of the negative of our function,
    #i.e. the maximum of our functino

    value = - value
    return value
end

function solveValueFunction(T, r, tol, minCons, numPointsA, Agrid)
    # will need to think of V1 and Agrid1
    ## ------------------------------------------------------------------------
    # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

    # Matrices to hold the policy and value functions
    V        = zeros(T+1, numPointsA)
    policyA1 = zeros(T,   numPointsA)
    policyC  = zeros(T,   numPointsA)

    ## ------------------------------------------------------------------------
    #Set the terminal value function to 0
    V[T + 1,:] = 0

    ## ------------------------------------------------------------------------
    # SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
    # BACKWARDS TO ZERO, ONE PERIOD AT A TIME

    for ixt=T:-1:1                                 # Loop from time T-1 to 1
        V1  = V[ixt + 1,:]                        #  Get tomorrow's value function
        Agrid1 = Agrid[ixt + 1, :]                # Get tomorrow's asset grid
        for ixA = 1:1:numPointsA                   # points on asset grid

            # Solve problem at grid points in assets
            # ---------------------------------------------------------
                # Information for optimisation
                A    = Agrid[ixt, ixA]            # assets today
                lbA1 = Agrid1[1]                  # lower bound: assets tomorrow
                ubA1 = (A - minCons)*(1+r)        # upper bound: assets tomorrow

                # define interpolation function itp
                knots_x = (Agrid1,)
                itp = interpolate(knots_x, V1, Gridded(Linear()))
                # @time itp[2.0]

                # Compute solution
                if (ubA1 - lbA1 < minCons)                         # if liquidity constrained
                    negV = objectivefunc(beta, r, itp, lbA1, A)
                    policyA1[ixt,ixA] = lbA1
                else
                    # define obj fcn
                    function obj(A1)
                        # println(A1)
                        return objectivefunc(beta, r, itp, A1, A)
                    end
                    Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          #Find the solution...
                    # println(Res)

                    policyA1[ixt,ixA] = Res.minimizer
                    negV = Res.minimum     # if interior solution
                end # if (ubA1 - lbA1 < minCons)
                #
                # # Store solution and its value
                policyC[ixt, ixA] = A - policyA1[ixt, ixA]/(1+r)
                V[ixt, ixA]       = -negV
        end #ixA

        # println(V[ixt,:])
        # println(policyC[ixt, :])
        println("Passed period $ixt of $T.")

    end #ixt

    return policyA1, policyC, V
end

function simNoUncer(T, r, Agrid, policyA1,V,startingA)
    # Initialise arrays that will hold the paths of income consumption, value
    # and assets

    # Arguments for output
    c = zeros(T, 1)            # consumption
    v = zeros(T, 1)            # value
    a = zeros(T + 1,1)         # this is the path at the start of each period, so we include the 'start' of death

    ## ------------------------------------------------------------------------
    # Obtain paths using the initial condition and the policy and value
    # functions
    #-------------------------------------------------------------------------#
    a[1, 1] = startingA
    for t = 1:1:T                     # loop through time periods for a particular individual
        # Original matlab code
        # v[t  , 1]   = interp1(Agrid[t, :],V[t, :],a[t, 1],interpMethod, 'extrap')
        knots_x = (Agrid[t, :],)
        itp = interpolate(knots_x, V[t, :], Gridded(Linear()))
        v[t  , 1]   = itp[ a[t, 1] ]

        # a[t+1, 1]   = interp1(Agrid[t, :], policyA1[t, :] ,a[t, 1],interpMethod, 'extrap')
        knots_x = (Agrid[t, :],)
        itp_assets  = interpolate(knots_x, policyA1[t, :], Gridded(Linear()))
        a[t+1, 1]   = itp_assets[ a[t, 1] ]

        c[t  , 1] = a[t, 1]  - (a[t+1, 1]/(1+r))
    end   #t
    return c, a, v
end

@time policyA1, policyC, V  = solveValueFunction(T, r, tol, minCons, numPointsA, Agrid)

cpath, apath, vpath = simNoUncer(T, r, Agrid, policyA1, V, startA);


# A = rand(20)
# A_x = collect(1.0:2.0:40.0)
# knots = (A_x,)
# @time itp = interpolate(knots, A, Gridded(Linear()))
# @time itp[2.0]
#
# print("done")

# Pkg.add("Plots")
# Pkg.add("PyPlot")
using Plots
pyplot() # Choose a backend

# ixt = 1
# x = Agrid[ixt, :]
# y = policyA1[ixt,:]
# plot(x,y)
#
# for ixt = 11:10:71
#     x = Agrid[ixt, :]
#     y = policyA1[ixt,:]
#     plot!(x,y)
# end
# gui()

plot(apath)
gui()
