
function solveValueFunction(minCons, Agrid)
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
                    negV = objectivefunc(itp, lbA1, A)
                    policyA1[ixt,ixA] = lbA1
                else
                    # define obj fcn
                    function obj(A1)
                        # println(A1)
                        return objectivefunc(itp, A1, A)
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

        println("Passed period $ixt of $T.")

    end #ixt

    return policyA1, policyC, V
end
