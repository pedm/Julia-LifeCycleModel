function solveEulerEquation(params::Dict{String,Float64}, Agrid, Ygrid, incTransitionMrx)
    # Define params
    minCons    = params["minCons"]
    r          = params["r"]

    # will need to think of V1 and Agrid1
    ## ------------------------------------------------------------------------
    # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

    # Matrices to hold the policy and value functions
    V        = zeros(T+1, numPointsA, numPointsY)
    policyA1 = zeros(T,   numPointsA, numPointsY)
    policyC  = zeros(T,   numPointsA, numPointsY)

    # Matrices to hold expected value and marginal utility functions
    EV  = zeros(T+1, numPointsA, numPointsY);
    dU  = zeros(T,   numPointsA, numPointsY);
    EdU = zeros(T,   numPointsA, numPointsY);

    ## ------------------------------------------------------------------------
    #Set the terminal value function and expected value function to 0
    V[T + 1, :, :] .= 0
    EV[T + 1, :, :] .= 0

    ## ------------------------------------------------------------------------
    # SOLVE THE HH PROBLEM FOR TIME T ONLY
    for ixt=T
        Agrid1 = Agrid[ixt + 1, :]                # The grid on assets tomorrow
        for ixA = 1:1:numPointsA                  # points on asset grid
            for ixY = 1:1:numPointsY              # points on income grid

                # Value of income and information for optimisation
                A    = Agrid[ixt, ixA]            # assets today
                Y    = Ygrid[ixt, ixY]            # income today

                # Store (imposed) solution and its value
                # policyC[ixt, ixA, ixY] =  A + Y # impose A1 = 0
                policyC[ixt, ixA, ixY] = max( A + Y, params["minCons"] ) # impose A1 = 0. ensure that there is never zero consumption, as that would result in Inf utility
                dU[ixt, ixA, ixY]      = getmargutility( params, policyC[ixt, ixA, ixY] )

                if saveValue_inEE
                    V[ixt, ixA, ixY]       = utility(params, policyC[T, ixA, ixY])
                end
            end #ixY

            # STEP 2. integrate out income today conditional on income
            # yesterday to get EV and EdU
            # --------------------------------------------------------
            realisedV  = V[ixt, ixA, :]
            realiseddU = dU[ixt, ixA, :]
            if ixt >= Tretire
                # when retired, cannot transition between income states (not that it matters when everyone gets zero each period)
                EV[ixt, ixA, :]  = realisedV
                EdU[ixt, ixA, :]  = realiseddU
            else
                # when working, can transition between income states
                for ixY = 1:1:numPointsY
                    EV[ixt, ixA, ixY]  = dot( incTransitionMrx[ixY,:], realisedV)
                    EdU[ixt, ixA, ixY] = dot( incTransitionMrx[ixY,:], realiseddU)
                end #ixY
            end

        end #ixA

        if ixt % 10 == 0
            println("Passed period $ixt of $T.")
        end
    end #ixt

    ## ------------------------------------------------------------------------
    # SOLVE RECURSIVELY BACKWARDS
    for ixt=(T-1):-1:1
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
                EdU1 = EdU[ixt + 1, :, ixY]       # relevant section of Edu matrix (in assets tomorrow)

                # Define interpolation function for expected (inverse) MU of assets tomorrow
                nodes = (Agrid1,)
                if linearise == 0
                    itp = interpolate(nodes, EdU1, Gridded(Linear()))
                elseif linearise == 1
                    linEdU1 = getinversemargutility(params,EdU1);
                    itp = interpolate(nodes, linEdU1, Gridded(Linear()));
                end

                signoflowerbound = sign(eulerforzero(params, itp, A, lbA1, Y))

                # Require lower bound to be <= 0
                if (signoflowerbound == 1.0) || (ubA1 - lbA1 < minCons)        # if liquidity constrained

                    policyA1[ixt,ixA,ixY] = lbA1

                else # if interior solution

                    # Require the upper bound and lower bound to be different signs
                    signofupperbound = sign( eulerforzero(params, itp, A, ubA1, Y) )
                    if (signoflowerbound*signofupperbound == 1.0)
                        # CRRA satisfies inada conditions so this should never happen
                        error("Sign of lower bound and upperbound are the same - no solution to Euler equation. Bug likely")
                    end

                    ## Root finding: custom bisection function
                    function ee(A1::Float64)
                        eulerforzero(params, itp, A, A1, Y)
                    end
                    policyA1[ixt, ixA, ixY] = bisection64_custom(ee, lbA1, ubA1, params["tol"])

                end # interior solution

                # Store solution and its value
                policyC[ixt, ixA, ixY] = A + Y - policyA1[ixt, ixA, ixY]/(1+r)
                dU[ixt, ixA, ixY]      = getmargutility( params, policyC[ixt, ixA, ixY] )

                if saveValue_inEE
                    # Save value function as well (not necessary for solution but might be wanted for e.g. welfare calcs)
                    EV1              = EV[ixt + 1,:, ixY] # relevant section of EV matrix (in assets tomorrow)
                    itpEV1              = interpolate(nodes, EV1, Gridded(Linear()))
                    V[ixt, ixA, ixY] = -objectivefunc(params, itpEV1, policyA1[ixt,ixA,ixY], A, Y)
                end
            end #ixY

            # STEP 2. integrate out income today conditional on income
            # yesterday to get EV and EdU
            # --------------------------------------------------------
            realisedV  = V[ixt, ixA, :]
            realiseddU = dU[ixt, ixA, :]
            for ixY = 1:1:numPointsY
                EV[ixt, ixA, ixY]  = dot( incTransitionMrx[ixY,:], realisedV)
                EdU[ixt, ixA, ixY] = dot( incTransitionMrx[ixY,:], realiseddU)
            end #ixY

        end #ixA

        #if ixt % 10 == 0
            println("Passed period $ixt of $T.")
        #end
    end #ixt

    return policyA1, policyC, V, EV, dU, EdU
end

