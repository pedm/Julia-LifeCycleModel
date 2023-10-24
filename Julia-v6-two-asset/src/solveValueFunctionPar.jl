function solveValueFunctionPar(params::Dict{String,Float64}, Agrid, Bgrid, Ygrid, incTransitionMrx)
    # Define params
    minCons     = params["minCons"]
    r           = params["r"]
    r_b         = params["r_b"]
    max_contrib = params["max_contrib"]
    
    # Simplifications (to define somewhere else automatically)
    R_a = 1.0+r
    R_b = 1.0+r_b
    R_b_over_R_a =  R_b / R_a
    R_a_over_R_b =  R_a / R_b
    numPointsAstar = numPointsA + 4 # depends on extra grid points added to allow Astar to go negative
    # numPointsAstar = numPointsA

    # will need to think of V1 and Agrid1
    ## ------------------------------------------------------------------------
    # GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS ZEROS

    # Matrices to hold the policy and value functions
    V     = zeros(T+1, numPointsA, numPointsB, numPointsY)
    policyA1    = zeros(T,   numPointsA, numPointsB, numPointsY)
    policyB1    = zeros(T,   numPointsA, numPointsB, numPointsY)
    policyC     = zeros(T,   numPointsA, numPointsB, numPointsY)

    # Objects to store things from the first stage "no adjust" solution (solving for liquid choice holding illiquid fixed)
    V_NA        = zeros(T+1, numPointsAstar, numPointsB, numPointsY) # Value fcn if not adjusting illiquid asset 
    policyA1_NA = zeros(T,   numPointsAstar, numPointsB, numPointsY) # Policy fcn A1 if not adjusting illiquid asset

    policyAdj     = zeros(T,   numPointsA, numPointsB, numPointsY)

    # Matrices to hold expected value and marginal utility functions
    EV  = zeros(T+1, numPointsA, numPointsB,numPointsY);

    ## ------------------------------------------------------------------------
    #Set the terminal value function and expected value function to 0
    V_NA[T + 1, :, :, :] .= 0 # value function conditional on not adjusting illiquid assets
    V[T + 1, :, :, :] .= 0       # value function with the adjustment
    EV[T + 1, :, :, :] .= 0      # expected value - the "post decision value function"

    # Allow astargrid to go negative
    # getGrid(minongrid, maxongrid, GridPoints, method)
    # astar_negs = getGrid(-max_contrib, 0.0, 10, "equalsteps")[1:end-1]

    ## ------------------------------------------------------------------------
    # SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
    # BACKWARDS TO ZERO, ONE PERIOD AT A TIME

    for ixt=T:-1:1                                # Loop from time T-1 to 1
        Agrid1 = Agrid[ixt + 1, :]                # The grid on assets tomorrow

        # TODO: need to define Astargrid over all potential middle-of-period liquid asset options -- then loop over that in stage 1
        # Then use this for interpolation in stage 2
        # Astargrid = Agrid[ixt, :] # for now, use this simple grid... later will need to allow it to go negative, if we have fixed adj cost
        # TODO: do i need to allow this to go negative?
        # Astargrid = [astar_negs; Agrid[ixt, :]] # for now, use this simple grid... later will need to allow it to go negative, if we have fixed adj cost
        # println(Astargrid)

        # STEP 1. Inner problem: solve for consumption assuming that the household does not adjust illiquid assets
        # ---------------------------------------------------------
        println("Step 1")
        Threads.@threads for ixB0 = 1:1:numPointsB      # points on illiq asset grid
        # for ixB0 = 1:1:numPointsB      # points on illiq asset grid
            Bgrid1 = Bgrid[ixt + 1, :]                # The grid on assets tomorrow
            B0     = Bgrid[ixt, ixB0]            # start-of-period illiq. assets

            for ixY = 1:numPointsY                      # points on income grid

                # define interpolation function itp - note does not depend on ixA0
                # nodes = (Agrid1,)
                # EV1  = EV[ixt + 1, :, ixB0, ixY]         # relevant section of EV matrix (in assets tomorrow)
                # itp = interpolate(nodes, EV1, Gridded(Linear()))
                # WARNING: currently I'm selecting B1 using ixB0... but the grid changes from one period to the next (!) so this doesn't actually maintain the same exact B1 = B0 !!!!
                # TODO: might need to do a 2D interpolation here!!!
                
                # Needed b/c Bgrid[ixt, ixB0] does not equal Bgrid[ixt+1, ixB0] -- might think later about some trick to do this
                nodes = (Agrid1,Bgrid1,)
                EV1  = EV[ixt + 1, :, :, ixY]         # relevant section of EV matrix (in assets tomorrow)
                itp = interpolate(nodes, EV1, Gridded(Linear()))

                Y0              = Ygrid[ixt, ixY]                                    # income today

                # Set astar_negs: IMPORTANT MUST BE SAME IN STAGE 1 AND STAGE 2 -- better to put this code elsewhere - though it does depend on Y0
                astar_min_uncon = (minCons-Y0-r_b*B0) / R_a # TODO: what is up with returns here?
                astar_min       = max(astar_min_uncon, -max_contrib)
                astar_negs      = getGrid(astar_min, 0.0, 5, "equalsteps")[1:end-1]
                Astargrid       = [astar_negs; Agrid[ixt, :]] 
                # Astargrid       = Agrid[ixt, :] 
                
                # maybe better to define Lgrid here and then loop over ixL0
                for ixA0 = 1:numPointsAstar               # points on liq asset grid

                    # Value of income and information for optimisation
                    A0    = Astargrid[ixA0]             # middle-of-period liquid assets
                    lbA1  = Agrid[ixt + 1, 1]            # lower bound: assets tomorrow
                    # ubA1 = (A0 + Y0 - minCons)*(1.0+r)  # OLD: upper bound: assets tomorrow
                    # ubA1   = R_a*A0 + Y0 - r_b*B0 - minCons  # upper bound: assets tomorrow
                    ubA1 = R_a*A0 + Y0 - minCons      # upper bound: assets tomorrow
                    # My note: pretty sure this should be + r_b*B0 !!

                    # WAIT... is it weird that you get your return on B0 now? When this is effectively a choice variable in the no-adjust problem?
                    # Maybe the return should already be actualized on what you were born with... not what you keep

                    # TODO: is ubA1 still accurate when we have the illiquid adjustment????

                    # Compute solution
                    if (ubA1 - lbA1 < minCons)        # if liquidity constrained
                        # WARNING: lots of places where ubA1 < 0... therefore even lbA1 not possible (due to -r_b*B0 above)
                        
                        # println([lbA1, ubA1])
                        # negV = objectivefunc(params, itp, lbA1, A0, Y0, B0)
                        negV = - ( utility(params, minCons) + params["beta"] * itp[lbA1, B0*R_b] )
                        policyA1_NA[ixt,ixA0,ixB0,ixY] = lbA1

                        # println("wait - does this part work?")
                    else
                        # Find the A1 that minimizes the objective function

                        # define obj fcn
                        function obj(A1::Float64)
                            return objectivefunc(params, itp, A1, A0, Y0, B0)
                        end
                        # try 
                        #     Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          # println(Res)
                        # catch
                        #     println("ERROR")
                        #     println([lbA1,ubA1])
                        #     Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          # println(Res)
                        # end

                        Res = optimize(obj,lbA1,ubA1, abs_tol = 1e-5)          # println(Res)
                        policyA1_NA[ixt,ixA0,ixB0,ixY] = Res.minimizer
                        negV = Res.minimum     # if interior solution
                    end # if (ubA1 - lbA1 < minCons)

                    # Store solution and its value
                    # policyC_NA[ixt, ixA0, ixB0, ixY] = A0 + Y0 - policyA1_NA[ixt, ixA0, ixB0, ixY]/(1.0+r)
                    # todo: policy fcn for deposit = 0 in this case 
                    V_NA[ixt, ixA0, ixB0, ixY]       = -negV
                end #ixA0
            end #ixY
        end # ixB0

        # The above gives V_NA(ixt, ixA0, ixB0, ixY) which assumes no adjustment of illiquid asset
        # Next step: find the optimal deposit d into the illiquid asset

        # STEP 2. Outer problem to find optimal illqiuid assets
        # Will have to interpoate V_NA over Lgrid and Bgrid 

        # Create Lgrid given ixB0 and ixY
        # Works because in the inner problem: a = l - c => can rearrange so that l = a + c       
        # Lgrid = Agrid + policyC_NA[ixt, :, ixB0, ixY]
        # But wait.... ixB0 can change?
        # I guess can loop over ixB0 and ixY, then for each create Lgrid, then interpolate V_NA on Lgrid and Bgrid

        Bgrid0 = Bgrid[ixt, :]

        println("Step 2")
        # STEP 2. Outer problem to find optimal illqiuid assets
        # ---------------------------------------------------------
        # Threads.@threads 
        for ixB0 = 1:1:numPointsB       # points on illiq asset grid
            for ixY = 1:numPointsY                       # points on income grid

                Y0              = Ygrid[ixt, ixY]  # income today
                B0              = Bgrid[ixt, ixB0] # start-of-period illiq. assets

                # Set astar_negs: IMPORTANT MUST BE SAME IN STAGE 1 AND STAGE 2 -- better to put this code elsewhere - though it does depend on Y0
                # BE VERY CAREFUL CHANGING THIS CODE SINCE IT SHOWS UP TWICE - AND FIX THIS LATER
                astar_min_uncon = (minCons-Y0-r_b*B0) / R_a 
                astar_min       = max(astar_min_uncon, -max_contrib)
                astar_negs      = getGrid(astar_min, 0.0, 5, "equalsteps")[1:end-1]
                Astargrid       = [astar_negs; Agrid[ixt, :]] 
                # Astargrid       = Agrid[ixt, :]

                # define interpolation function over V_NA (as a function of middle-of-period liquid assets = a_star in Seb Graves language)
                nodes = (Astargrid, Bgrid0,)
                V_NA_slice  = V_NA[ixt, :, :, ixY]         # relevant section of EV matrix (in assets tomorrow)
                itp_V_NA = interpolate(nodes, V_NA_slice, Gridded(Linear()))

                # TODO: try without this!
                itp_V_NA = extrapolate(itp_V_NA, Line()) # Added because sometimes astar can go above astar grid... but maybe we want to think about this more later...
                
                for ixA0 = 1:numPointsA               # points on liq asset grid
                    A0     = Agrid[ixt, ixA0]         # start-of-period liquid assets

                    # Define a_star as middle-of-period liquid assets (i.e. what you have after adjusting your illiquid assets)
                    function a_star(B1) # middle-of-period liquid assets
                        return A0 + R_b_over_R_a*(B0 - B1) - ( transaction_costs(ixt, B1, B0) / R_a )
                    end

                    # define obj fcn (Seb Graves' Vtilde)
                    function Vtilde_obj(B1::Float64)
                        a_middle_of_period = a_star(B1)

                        # if a_middle_of_period < 0
                        #     println(B1)
                        #     println(B0)
                        #     println(a_middle_of_period)
                        #     error("too low")
                        # end
                        
                        try
                            return - itp_V_NA(a_middle_of_period, B1)
                        catch 
                            println("found issue")
                            println(B1)
                            println(a_middle_of_period)
                            println(nodes)
                            println(V_NA_slice)
                            return - itp_V_NA(a_middle_of_period, B1)
                        end
                    end

                    # Next step: choose B1 to get optimal Vtilde(B1) subject to constraints

                    # Choose upper and lower bound for illiquid assets tomorrow
                    lbB1     = Bgrid[ixt + 1, 1]            # NOTE: later might want to only allow a certain amount of withdrawal
                    ubB1_max = B0*(1.0 + params["r_b"]) + min(params["max_contrib"], (A0 + Y0 - minCons)*(1.0+r)) # current illiquid assets plus max contribution or total liq resources, whichever is smaller
                    # ubB1_max = B0*(1.0 + params["r_b"]) + (A0 + Y0 - minCons)*(1.0+r) # current illiquid assets plus total liq resources
                    # NOTE: or better to put the constraint afterwards? like, first do the uncon problem for B1, then impose constraints? dunno 

                    if (ubB1_max - lbB1 < minCons)        # issue - look into this!
                        println(ubB1_max)
                        println(lbB1)
                        error("some issue where ubB1_max < lbB1... why?")
                    end

                    # println( (lbB1, B0, ubB1_max) )
                    # println( (a_star(lbB1), a_star(B0), a_star(ubB1_max)) )
                    # println("")

                    # Sometimes the discreteness of the adjustment cost means that a_star(lbB1) < 0
                    # In that case, don't allow them to withdraw from their retirement account
                    if a_star(lbB1) < 0
                        lbB1 = B0 
                    end


                    # if a_star(lbB1) * a_star(ubB1_max) >= 0
                    #     println("ISSUE")
                    #     println(lbB1)
                    #     println( a_star(lbB1) )
                    #     println(ubB1_max)
                    #     println( a_star(ubB1_max) )

                    #     B1opts = collect(lbB1:.01:ubB1_max)
                    #     println( B1opts )
                    #     println( [a_star(b1) for b1 in B1opts] )
                    # end

                    # Find B1_overbar using analytical solution 
                    # Designed to ensure that astar >= astar_min which in turn ensures that A1_NA >= 0 and C_NA > 0 
                    B1_overbar = R_a_over_R_b * (A0 - astar_min) + B0 

                    # Find ubB1: the max illiquid assets such that a_star >= zero
                    # ubB1 = find_zero(a_star, (lbB1, ubB1_max), Bisection())
                    ubB1 = min(ubB1_max, B1_overbar)
                    # Note that a_star is decreasing in B1, if above B0, as putting more into retirement means less in liquid
                    # That said, a_star may decrease just below B0 due to the adjustment cost

                    # TODO: might want to have some cap on retirement contributions each period, no? 
                    # TODO: there is probably a much cleaner way to do this analytically, no?

                    # Maybe better: a_star of zero -> back out B1 -> make that the max -> then see if that's larger than max contribution / max assets
                    # But wait... why do we need to ensure that a_star is zero? Should we not allow a_star to go negative??

                    # TODO: GOT RID OF THIS, BUT SHOULD LOOK INTO IT
                    # if a_star(ubB1) + 0.000000001 < minimum(Astargrid)
                    #     println(ubB1)
                    #     println(a_star(ubB1))
                    #     println(astar_min)

                    #     println([A0, B0, Y0])
                    #     println([ubB1_max, B1_overbar])

                    #     println(Astargrid)

                    #     println("this will be an issue")
                    #     error("ugh")
                    # end

                    if isapprox(ubB1, lbB1)
                        # IF UNABLE TO ADJUST: just set equal to lbB1
                        policyB1[ixt,ixA0,ixB0,ixY]       = lbB1
                        V[ixt, ixA0, ixB0, ixY]           = - V_NA[ixt, ixA0, ixB0, ixY]
                        policyA1[ixt,ixA0,ixB0,ixY]       = policyA1_NA[ixt, ixA0, ixB0, ixY]
                        # Wait.... on the RHS should that be ixB0 = 1??
                        # Also wait... in this case, shouldn't we feed in lbB1 into the interpolated functions?
                        
                        println("ubB1 = lbB1 weird!!! ")
                    elseif (ubB1 - lbB1 < minCons)        # issue - look into this!
                        # IF ubB1 is less than lbB1... werid... 

                        println(ixt)
                        println(A0)
                        println(B0) 
                        println(Y0)

                        println(ubB1)
                        println(lbB1)
                        println(a_star(ubB1))
                        println(a_star(lbB1))
                        println(minCons)
                        error("some issue where ubB1 < lbB1... why?")
                    else
                        # IF ABLE TO ADJUST
                        # NOTE: here we always assume interior solution... not sure if that's true
                        # println([A0, B0, Y0, lbB1, ubB1, a_star(ubB1), maximum(Astargrid), a_star(ubB1) > maximum(Astargrid)])
                        # println([A0, B0, Y0, lbB1, ubB1 ])

                        # try
                        #     Res = optimize(Vtilde_obj, lbB1, ubB1, abs_tol = 1e-5)          # println(Res)
                        # catch 
                        #     println("")
                        #     println("CATCHING ERROR")
                        #     println([lbB1, ubB1])
                        #     b1s = [collect(lbB1:0.1:ubB1); ubB1]
                        #     a_stars = [a_star(b1) for b1 in b1s]
                        #     println(b1s)
                        #     println(a_stars)
                        #     plot(b1s, a_stars)
                        #     gui()

                        #     Res = optimize(Vtilde_obj, lbB1, ubB1, abs_tol = 1e-5)          # println(Res)
                        # end

                        Res = optimize(Vtilde_obj, lbB1, ubB1, abs_tol = 1e-5)          # println(Res)
                        policyB1[ixt,ixA0,ixB0,ixY]     = Res.minimizer
                        V[ixt, ixA0, ixB0, ixY]         = - Res.minimum     
                        policyAdj[ixt, ixA0, ixB0, ixY] = 1.0

                        # Interpolate over policy function for liquid assets assuming no adjustment of B1
                        # define interpolation function over V_NA (as a function of middle-of-period liquid assets = a_star in Seb Graves language)
                        policyA1_NA_slice  = policyA1_NA[ixt, :, :, ixY]         # relevant section of EV matrix (in assets tomorrow)
                        itp_policyA1_NA = interpolate(nodes, policyA1_NA_slice, Gridded(Linear()))
                        itp_policyA1_NA = extrapolate(itp_policyA1_NA, Line()) # Added because sometimes astar can go below astar grid... but maybe we want to think about this more later...

                        # TODO: alternative here....
                        # Rather than interpolate over V_NA... could just use the policy function for A1 given astar and B1 -> then use that to get consumption -> then based on that get u(c)+beta*EV(.)

                        # Set A1 policy fcn based on chosen B1 and the mapping from B1->A*->A1
                        try
                            policyA1[ixt,ixA0,ixB0,ixY] = itp_policyA1_NA( a_star( policyB1[ixt,ixA0,ixB0,ixY] ), policyB1[ixt,ixA0,ixB0,ixY] ) # See bottom of page 5 from Seb Graves notes
                        catch
                            println("Issue found:")
                            println([ixt,ixA0,ixB0,ixY])
                            println([ixt, A0, B0, Y0])
                            println([lbB1, ubB1, Res.minimizer])
                            println([a_star(lbB1), a_star(ubB1), a_star(Res.minimizer)])
                            println("Im confused:")
                            println( A0 + R_b_over_R_a*(B0 - Res.minimizer) )
                            println( R_b_over_R_a*(B0 - Res.minimizer)  )
                            println(Astargrid)
                            policyA1[ixt,ixA0,ixB0,ixY] = itp_policyA1_NA( a_star( policyB1[ixt,ixA0,ixB0,ixY] ), policyB1[ixt,ixA0,ixB0,ixY] ) # See bottom of page 5 from Seb Graves notes

                            # Going over astar grid -- this surprises me. i would have thought that if B1 = 0, then Astar = A0 should hold. But it doesnt!????
                            # Oh it's because B0 = 0.01 and B1 = 0.0 -> thus we're over the max Astar grid. Shit. What to do about this? 
                            # Hmmm why did it choose to extract in that case? Something to do with extrapolation above? Not sure. 
                            # I guess the lbB1 is always just 0... as opposed to something that ensures we dont end up above the Agrid. 
                            # Hmm conceptually, what should we do here? 
                        end
                    end
                end #ixA0
            end #ixY
        end # ixB0

        # Notes from Jeppe: 
        # m_t+1 = R * a_t + y_t+1

        # Or easier for VFI:
        # interpolate V_NA over Agrid and Bgrid
        # Then use the fact that a = beg period liq assets + whatever you extract - deposits

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
        # if ixt % 10 == 0
            println("Passed period $ixt of $T.")
        # end
    end #ixt

    return policyA1, policyB1, policyC, V, EV, V_NA, policyAdj
end
