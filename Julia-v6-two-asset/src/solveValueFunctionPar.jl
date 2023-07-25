function solveValueFunctionPar(params::Dict{String,Float64}, Agrid, Bgrid, Ygrid, incTransitionMrx)
    # Define params
    minCons    = params["minCons"]
    r          = params["r"]

    # will need to think of V1 and Agrid1
    ## ------------------------------------------------------------------------
    # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS ZEROS

    # Matrices to hold the policy and value functions
    V_inner  = zeros(T+1, numPointsA, numPointsB, numPointsY)
    V        = zeros(T+1, numPointsA, numPointsB, numPointsY)
    policyA1_inner = zeros(T,   numPointsA, numPointsB, numPointsY)
    policyA1 = zeros(T,   numPointsA, numPointsB, numPointsY)
    policyC_inner  = zeros(T,   numPointsA, numPointsB, numPointsY)
    policyC  = zeros(T,   numPointsA, numPointsB, numPointsY)

    # Matrices to hold expected value and marginal utility functions
    EV  = zeros(T+1, numPointsA, numPointsB,numPointsY);

    # Liquid assets just before consumption 
    Lgrid = zeros(size(Agrid))

    ## ------------------------------------------------------------------------
    #Set the terminal value function and expected value function to 0
    V_inner[T + 1, :, :, :] .= 0 # value function conditional on not adjusting illiquid assets
    V[T + 1, :, :, :] .= 0       # value function with the adjustment
    EV[T + 1, :, :, :] .= 0      # expected value - the "post decision value function"

    ## ------------------------------------------------------------------------
    # SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
    # BACKWARDS TO ZERO, ONE PERIOD AT A TIME

    for ixt=T:-1:1                                # Loop from time T-1 to 1
        Agrid1 = Agrid[ixt + 1, :]                # The grid on assets tomorrow

        # STEP 1. Inner problem: solve for consumption assuming that the household does not adjust illiquid assets
        # ---------------------------------------------------------
        Threads.@threads for ixB = 1:1:numPointsB       # points on illiq asset grid

            for ixY = 1:numPointsY                      # points on income grid

                # define interpolation function itp - note does not depend on ixA
                nodes = (Agrid1,)
                EV1  = EV[ixt + 1, :, ixB, ixY]         # relevant section of EV matrix (in assets tomorrow)
                itp = interpolate(nodes, EV1, Gridded(Linear()))
                
                for ixA = 1:numPointsA               # points on liq asset grid

                    # Value of income and information for optimisation
                    A    = Agrid[ixt, ixA]            # assets today
                    Y    = Ygrid[ixt, ixY]            # income today
                    lbA1 = Agrid[ixt + 1, 1]          # lower bound: assets tomorrow
                    ubA1 = (A + Y - minCons)*(1.0+r)  # upper bound: assets tomorrow

                    # Compute solution
                    if (ubA1 - lbA1 < minCons)        # if liquidity constrained
                        negV = objectivefunc(params, itp, lbA1, A, Y)
                        policyA1_inner[ixt,ixA,ixB,ixY] = lbA1
                    else
                        # Find the A1 that minimizes the objective function

                        # define obj fcn
                        function obj(A1::Float64)
                            return objectivefunc(params, itp, A1, A, Y)
                        end
                        Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          # println(Res)
                        policyA1_inner[ixt,ixA,ixB,ixY] = Res.minimizer
                        negV = Res.minimum     # if interior solution
                    end # if (ubA1 - lbA1 < minCons)

                    # Store solution and its value
                    policyC_inner[ixt, ixA, ixB, ixY] = A + Y - policyA1_inner[ixt, ixA, ixB, ixY]/(1.0+r)
                    # todo: policy fcn for deposit = 0 in this case 
                    V_inner[ixt, ixA, ixB, ixY]       = -negV
                end #ixA
            end #ixY
        end # ixB

        # Basically the above gives V_t(ixA, ixB, ixY) assuming no adjustment of illiquid asset
        # Next step: find the optimal deposit d into the illiquid asset
        # This will reduce 

        # STEP 2. Outer problem to find optimal deposit function
        # Will have to interpoate V_inner over Lgrid and Bgrid 
        
        # Create Lgrid given ixB and ixY
        # Works because in the inner problem: a = l - c => can rearrange so that l = a + c       
        # Lgrid = Agrid + policyC_inner[ixt, :, ixB, ixY]
        # But wait.... ixB can change?

        # Or easier for VFI:
        # interpolate V_inner over Agrid and Bgrid
        # Then use the fact that a = beg period liq assets + whatever you extract - deposits

        V[ixt, :, :, :] = V_inner[ixt, :, :, :]
        policyC[ixt, :, :, :] = policyC_inner[ixt, :, :, :]
        policyA1[ixt, :, :, :] = policyA1_inner[ixt, :, :, :]

        # STEP 3. integrate out income today conditional on income yesterday to get EV and EdU
        # "Post Decision Value Function" in the words of Jeppe 
        # --------------------------------------------------------
        Threads.@threads for ixA = 1:1:numPointsA                  # points on asset grid
            for ixB = 1:numPointsB
                realisedV = V[ixt, ixA, ixB, :]
                for ixY = 1:1:numPointsY
                    EV[ixt, ixA, ixB, ixY]  = dot(incTransitionMrx[ixY,:], realisedV)
                end #ixY
            end #ixB
        end #ixA

        # println(V[ixt,:])
        # println(policyC[ixt, :])
        if ixt % 10 == 0
            println("Passed period $ixt of $T.")
        end
    end #ixt

    return policyA1, policyC, V, EV
end