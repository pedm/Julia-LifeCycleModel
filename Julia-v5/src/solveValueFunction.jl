function solveValueFunction(Agrid, Ygrid, incTransitionMrx)
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
        Agrid1 = Agrid[ixt + 1, :]                # The grid on assets tomorrow

        for ixA = 1:1:numPointsA                  # points on asset grid

            # STEP 1. solve problem at grid points in assets and income
            # ---------------------------------------------------------
            for ixY = 1:1:numPointsY               # points on income grid

                # Value of income and information for optimisation
                A    = Agrid[ixt, ixA]            # assets today
                Y    = Ygrid[ixt, ixY]            # income today
                lbA1 = Agrid[ixt + 1, 1]          # lower bound: assets tomorrow
                ubA1 = (A + Y - minCons)*(1+r)    # upper bound: assets tomorrow
                EV1  = EV[ixt + 1,:, ixY]         # relevant section of EV matrix (in assets tomorrow)

                # define interpolation function itp
                knots_x = (Agrid1,)
                itp = interpolate(knots_x, EV1, Gridded(Linear()))

                # TODO: why so much slower with float rather than int?
                # @time itp[2.0]
                # @time itp[2.0]
                #
                # @time itp[2]
                # @time itp[2]
                # hgdgfd
                # Compute solution
                if (ubA1 - lbA1 < minCons)        # if liquidity constrained
                    negV = objectivefunc(itp, lbA1, A, Y)
                    policyA1[ixt,ixA,ixY] = lbA1
                else
                    # Find the A1 that minimizes the objective function

                    # define obj fcn
                    function obj(A1)
                        return objectivefunc(itp, A1, A, Y)
                    end
                    Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          # println(Res)
                    policyA1[ixt,ixA,ixY] = Res.minimizer
                    negV = Res.minimum     # if interior solution
                end # if (ubA1 - lbA1 < minCons)

                # Store solution and its value
                policyC[ixt, ixA, ixY] = A + Y - policyA1[ixt, ixA, ixY]/(1+r)
                V[ixt, ixA, ixY]       = -negV
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
