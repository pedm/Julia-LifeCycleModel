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
    V[T + 1, :, :] = 0
    EV[T + 1, :, :] = 0

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
            
function  getmargutility(cons)
   if cons<=0
                     error("Consumption is <=0")
      
   elseif gamma==1          
                    margut=1./cons
                    
   else
                    margut = (cons).^(-gamma)      
   end
end               
            
function getinversemargutility(margut)
   if gamma == 1
                 invmargut = 1/margut;
   elseif
                 invmargut = margut^(-1/gamma);
   end
end    

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
                # EV1  = EV[ixt + 1,:, ixY]       # relevant section of EV matrix (in assets tomorrow)
                EdU1 = EdU[ixt + 1, :, ixY]        # relevant section of Edu matrix (in assets tomorrow)

                # Define interpolation function itp
                 if linearise == 0
                 du1AtA1 = interp1(Agrid1,Edu1,A1, interpMethod, 'extrap');
                    knots_x = (Agrid1,)
                    EdU1_at_A1 = interpolate(knots_x, EdU1, Gridded(Linear()))
                 elseif linearise == 1
                     linEdU1 = getinversemargutility(Edu1);
                     invDu1atA1 = interp1(Agrid1,linEdU1,A1, interpMethod, 'extrap');
                     du1AtA1 = getmargutility(invDu1atA1);
                 end

                signoflowerbound = sign(eulerforzero(params, EdU1_at_A1, A, lbA1, Y))

                # Require lower bound to be <= 0
                if (signoflowerbound == 1.0) || (ubA1 - lbA1 < minCons)        # if liquidity constrained
                    # TODO: why would a positive signoflowerbound indicate liq
                    # constrained? (I understand the second part of this statement,
                    # but not the first)

                    # negV = objectivefunc(params, itp, lbA1, A, Y)
                    policyA1[ixt,ixA,ixY] = lbA1

                else # if interior solution

                    # Require the upper bound and lower bound to be different signs
                    signofupperbound = sign( eulerforzero(params, EdU1_at_A1, A, ubA1, Y) )
                    if (signoflowerbound*signofupperbound == 1.0)
                        # TODO why would this indicate a bug?
                        error("Sign of lower bound and upperbound are the same - no solution to Euler equation. Bug likely")
                    end

                    ############################################################
                    ## Root finding option 1: use the Roots package.
                    ############################################################
                    # Problem with the Roots package: Bisection does not allow you to set the tolerance
                    # function ee(A1::Float64)
                    #     eulerforzero(params, EdU1_at_A1, A, A1, Y)
                    # end
                    # policyA1[ixt, ixA, ixY] = find_zero( ee, [lbA1, ubA1], Bisection())

                    ############################################################
                    ## Root finding option 2: custom bisection function
                    ############################################################
                    # New version of bisection -> now i can control tolerance
                    function ee(A1::Float64)
                        eulerforzero(params, EdU1_at_A1, A, A1, Y)
                    end
                    policyA1[ixt, ixA, ixY] = bisection64_custom(ee, lbA1, ubA1, 0.001)

                    ############################################################
                    ## Root finding option 3: custom bisection function, with argument passed directly
                    ############################################################
                    # Same as bisection64_custom, but now i pass args directly
                    # NOTE: seems it doesnt make any difference
                    # policyA1[ixt, ixA, ixY] = bisection64_with_args(eulerforzero_rearrange, lbA1, ubA1, 0.001, (params, EdU1_at_A1, A, Y))

                end # if (ubA1 - lbA1 < minCons)

                # Store solution and its value
                policyC[ixt, ixA, ixY] = A + Y - policyA1[ixt, ixA, ixY]/(1+r)
                dU[ixt, ixA, ixY]      = getmargutility( params, policyC[ixt, ixA, ixY] )

                if saveValue_inEE
                    # Save value function as well (not strictly necessary if i do simulation using policy functions)
                    EV1              = EV[ixt + 1,:, ixY] # relevant section of EV matrix (in assets tomorrow)
                    knots_x          = (Agrid1,)
                    itp              = interpolate(knots_x, EV1, Gridded(Linear()))
                    V[ixt, ixA, ixY] = utility(params, policyC[ixt, ixA, ixY]) + params["beta"] * itp[ policyA1[ixt,ixA,ixY] ]
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

        if ixt % 10 == 0
            println("Passed period $ixt of $T.")
        end
    end #ixt

    return policyA1, policyC, V, EV, dU, EdU
end

# I pass in the args to f for speed, as suggested here:
# https://mmas.github.io/bisection-method-julia
function bisection64_with_args(f, a::Float64, b::Float64, tol::Float64, args=()::Tuple )

    if a > b
        b,a = a, b
    end

    m = _middle(a,b)
    fa, fb = sign(f(a, args...)), sign(f(b, args...))

    fa * fb > 0 && throw(ArgumentError(bracketing_error))
    (iszero(fa) || isnan(fa) || isinf(fa)) && return a
    (iszero(fb) || isnan(fb) || isinf(fb)) && return b

    while a < m < b
        f_val = f(m, args...)
        fm = sign(f_val)
        # println("m = $m and f_val = $f_val")

        if (abs(f_val) < tol) || iszero(fm) || isnan(fm) || isinf(fm)
            return m
        elseif fa * fm < 0
            b,fb=m,fm
        else
            a,fa=m,fm
        end
        m = _middle(a,b)
    end
    return m
end

function eulerforzero_rearrange(A1::Float64, params::Dict{String,Float64}, EdU1_at_A1, A0::Float64, Y::Float64)
    eulerforzero(params, EdU1_at_A1, A0, A1, Y)
end
