function solveValueFunction(params::Dict{String,Float64}, Agrid, Ygrid, incTransitionMrx)


    # will need to think of V1 and Agrid1
    ## ------------------------------------------------------------------------
    # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

    # Matrices to hold the policy and value functions
    V        = zeros(T+1, numPointsA, numPointsY)
    policyA1 = zeros(T,   numPointsA, numPointsY)
    policyC  = zeros(T,   numPointsA, numPointsY)

    # Matrices to hold expected value and marginal utility functions
    EV  = zeros(T+1, numPointsA, numPointsY);
    # EdU = zeros(T,   numPointsA, numPointsY);

    ## ------------------------------------------------------------------------
    #Set the terminal value function and expected value function to 0
    V[T + 1, :, :] = 0
    EV[T + 1, :, :] = 0

    ## ------------------------------------------------------------------------
    # SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
    # BACKWARDS TO ZERO, ONE PERIOD AT A TIME

    for ixt=T:-1:1                                # Loop from time T-1 to 1

        for ixA = 1:1:numPointsA                  # points on asset grid

            # STEP 1. solve problem at grid points in assets and income
            # ---------------------------------------------------------
            for ixY = 1:1:numPointsY               # points on income grid
                # relevant section of EV matrix (in assets tomorrow)
                # NOTE: that EV1 will only depend on the persistent component of income
                EV1                                                 = EV[ixt + 1,:, ixY]

                # Assuming no transitory income shocks
                # negV, policyA1[ixt,ixA,ixY], policyC[ixt, ixA, ixY] = Value(params, Agrid, Ygrid, EV1, ixt, ixA, ixY, 0.0)

                # Allowing for transitory income shocks
                function Value_times_pdf(Y_trans::Float64)
                    return Value_return1(params, Agrid, Ygrid, EV1, ixt, ixA, ixY, Y_trans) * normpdf(-0.5 * params["sigma_trans"], params["sigma_trans"], Y_trans)
                end

                # NOTE: problem, this integration method will not work because we want to obtain policy functions too
                negV = quadgk(Value_times_pdf, -normBnd * params["sigma_trans"], normBnd * params["sigma_trans"]; reltol = 0.1, abstol=0.1)[1]
                # QUESTION: is it acceptable to cutoff at 0.0 rather than -normBnd * σ ???

                # TODO: will need to produce policy functions!!!!

                V[ixt, ixA, ixY]                                    = -negV
            end #ixY

            # STEP 2. integrate out income today conditional on income
            # yesterday to get EV and EdU
            # --------------------------------------------------------
            realisedV = V[ixt, ixA, :]
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
    Y    = Ygrid[ixt, ixY] + Y_trans           # income today
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
