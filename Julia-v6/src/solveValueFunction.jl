function solveValueFunction(params::Dict{String,Float64}, Agrid, Ygrid, incTransitionMrx)


    # will need to think of V1 and Agrid1
    ## ------------------------------------------------------------------------
    # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

    # Matrices to hold the policy and value functions
    V        = zeros(T+1, numPointsA, numPointsY, numPointsYTrans)
    policyA1 = zeros(T,   numPointsA, numPointsY, numPointsYTrans)
    policyC  = zeros(T,   numPointsA, numPointsY, numPointsYTrans)

    # Setup nodes and weights for transitory income shocks
    μtransshocks = -0.5 * params["sigma_trans"] # this comes from Zeldes 1989 and Kovacs 2015. The expected value of a log normal variable with mean μ and variance σ2 is given by exp(μ + σ2/2)
    Ytrans_grid, Ytrans_weights = gausshermite_normal_distribution(numPointsYTrans, μtransshocks, params["sigma_trans"] )
    Ytrans_grid = exp.(Ytrans_grid)

    # Matrices to hold expected value and marginal utility functions
    EV  = zeros(T+1, numPointsA, numPointsY);
    # EdU = zeros(T,   numPointsA, numPointsY);

    ## ------------------------------------------------------------------------
    #Set the terminal value function and expected value function to 0
    V[T + 1, :, :, :] = 0
    EV[T + 1, :, :] = 0

    ## ------------------------------------------------------------------------
    # SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
    # BACKWARDS TO ZERO, ONE PERIOD AT A TIME

    for ixt=T:-1:1                                # Loop from time T-1 to 1

        for ixA = 1:1:numPointsA                  # points on asset grid

            # STEP 1. solve problem at grid points in assets and income
            # ---------------------------------------------------------
            for ixY = 1:1:numPointsY               # points on income grid
                for ixYtrans = 1:numPointsYTrans
                    # relevant section of EV matrix (in assets tomorrow)
                    # NOTE: that EV1 will only depend on the persistent component of income
                    EV1                                                 = EV[ixt + 1,:, ixY]

                    # Assuming no transitory income shocks
                    # negV, policyA1[ixt,ixA,ixY], policyC[ixt, ixA, ixY] = Value(params, Agrid, Ygrid, EV1, ixt, ixA, ixY, 0.0)

                    # Allowing for transitory income shocks
                    # function Value_times_pdf(Y_trans::Float64)
                    #     return Value_return1(params, Agrid, Ygrid, EV1, ixt, ixA, ixY, Y_trans) * normpdf(μtransshocks, params["sigma_trans"], Y_trans)
                    # end
                    # negV = quadgk(Value_times_pdf, -normBnd * params["sigma_trans"], normBnd * params["sigma_trans"]; reltol = 0.1, abstol=0.1)[1]
                    # NOTE: problem, this integration method will not work because we want to obtain policy functions too
                    # TODO: will need to produce policy functions!!!!

                    Ytrans = Ytrans_grid[ixYtrans]

                    negV, policyA1[ixt,ixA,ixY, ixYtrans], policyC[ixt, ixA, ixY, ixYtrans] = Value(params, Agrid, Ygrid, EV1, ixt, ixA, ixY, Ytrans)
                    V[ixt, ixA, ixY, ixYtrans]                                              = -negV
                end #ixY_trans
            end #ixY

            # STEP 2. integrate out income today conditional on income
            # yesterday to get EV and EdU
            # --------------------------------------------------------
            # Integrate over transitory income shocks
            realisedV = V[ixt, ixA, :, :] * Ytrans_weights
            # Integrate over persistent income shocks conditional on previous ixY
            for ixY = 1:1:numPointsY
                EV[ixt, ixA, ixY]  = dot( incTransitionMrx[ixY,:], realisedV)
            end #ixY

        end #ixA

        # println(V[ixt,:])
        # println(policyC[ixt, :])
        if ixt % 10 == 0
            println("Passed period $ixt of $T.")
        end
    end #ixt

    return policyA1, policyC, V, EV
end


function Value(params::Dict{String,Float64}, Agrid, Ygrid, EV1, ixt, ixA, ixY, Y_trans)

    # Define params and model objects
    minCons    = params["minCons"]
    r          = params["r"]
    Agrid1     = Agrid[ixt + 1, :]                # The grid on assets tomorrow

    # Value of income and information for optimisation
    A    = Agrid[ixt, ixA]                     # assets today
    Y    = Ygrid[ixt, ixY] * Y_trans           # income today

    # Government imposed income floor
    if Y < minCons
        Y = minCons
    end

    lbA1 = Agrid[ixt + 1, 1]                   # lower bound: assets tomorrow
    ubA1 = (A + Y - minCons)*(1+r)             # upper bound: assets tomorrow

    # define interpolation function itp
    knots_x = (Agrid1,)
    itp = interpolate(knots_x, EV1, Gridded(Linear()))

    # TODO: why so much slower with float rather than int?
    # @time itp[2.0]
    # @time itp[2.0]
    #
    # @time itp[2]
    # @time itp[2]

    # Compute solution
    if (ubA1 - lbA1 < minCons)        # if liquidity constrained
        negV = objectivefunc(params, itp, lbA1, A, Y)
        policyA1 = lbA1
    else
        # Find the A1 that minimizes the objective function

        # define obj fcn
        function obj(A1::Float64)
            return objectivefunc(params, itp, A1, A, Y)
        end
        Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          # println(Res)
        policyA1 = Res.minimizer
        negV = Res.minimum     # if interior solution
    end # if (ubA1 - lbA1 < minCons)

    # Store solution and its value
    policyC = A + Y - policyA1/(1+r)

    return negV, policyA1, policyC
end

function Value_return1(params::Dict{String,Float64}, Agrid, Ygrid, EV1, ixt, ixA, ixY, Y_trans)

    # Define params and model objects
    minCons    = params["minCons"]
    r          = params["r"]
    Agrid1     = Agrid[ixt + 1, :]                # The grid on assets tomorrow

    # Value of income and information for optimisation
    A    = Agrid[ixt, ixA]                     # assets today
    Y    = Ygrid[ixt, ixY] + Y_trans           # income today
    # Government imposed income floor
    if Y < minCons
        Y = minCons
    end
    lbA1 = Agrid[ixt + 1, 1]                   # lower bound: assets tomorrow
    ubA1 = (A + Y - minCons)*(1+r)             # upper bound: assets tomorrow

    # define interpolation function itp
    knots_x = (Agrid1,)
    itp = interpolate(knots_x, EV1, Gridded(Linear()))

    # TODO: why so much slower with float rather than int?
    # @time itp[2.0]
    # @time itp[2.0]
    #
    # @time itp[2]
    # @time itp[2]

    # Compute solution
    if (ubA1 - lbA1 < minCons)        # if liquidity constrained
        negV = objectivefunc(params, itp, lbA1, A, Y)
        policyA1 = lbA1
    else
        # Find the A1 that minimizes the objective function

        # define obj fcn
        function obj(A1::Float64)
            return objectivefunc(params, itp, A1, A, Y)
        end
        Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          # println(Res)
        policyA1 = Res.minimizer
        negV = Res.minimum     # if interior solution
    end # if (ubA1 - lbA1 < minCons)

    # Store solution and its value
    policyC = A + Y - policyA1/(1+r)

    return negV
end

# Gauss Hermite Quadrature - meant for integrating from -Inf to Inf
# where the inside of the integral is something of the form e^(-x2) * f(x)
# https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature
# https://jblevins.org/notes/quadrature
# Note that we have to change the variable slightly to get the normal distribution into that form
# Say we want to integrate ∫ p(x) V(x) dx where p(x) is the normal pdf and V(x) is the value function evaluated at x

function gausshermite_normal_distribution(n, μ, σ)
    nodes, weights = gausshermite( n ) # from package FastGaussQuadrature
    nodes_mod      = sqrt(2.0) * σ .* nodes + μ
    weights_mod    = weights ./ sqrt(π)
    return nodes_mod, weights_mod
end

# NOTE that the nodes will sum up to 1. This is the same as integrating over the normal pdf from -Inf to Inf
# nodes, weights = gausshermite_normal_distribution(9, 0.0, sqrt(0.05))
