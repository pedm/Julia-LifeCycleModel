

function simWithUncer(model, Agrid, Bgrid, Ygrid, policyA1, policyB1, EV)
    # (policyA1,EV,startA)

    # This function takes the policy functions and value functions, along with
    # starting assets and returns simulated paths of income, consumption,
    # assets and value

    params     = model["params"]
    ints       = model["ints"]
    # T          = ints["T"]
    # Tretire    = ints["Tretire"]
    # numPointsY = ints["numPointsY"]
    # numPointsA = ints["numPointsA"]
    # numPointsB = ints["numPointsB"]
    numSims    = ints["numSims"]

    ## ------------------------------------------------------------------------
    # Define parameters

    rho     = params["rho"]
    mu      = params["mu"]
    sigma   = params["sigma"]
    r       = params["r"]
    startA  = params["startA"]

    ## ------------------------------------------------------------------------
    # Initialise arrays that will hold the paths of income consumption, value
    # and assets

    # Arguments for output
    y = zeros(T, numSims);            # income
    c = zeros(T, numSims);            # consumption
    ew= zeros(T, numSims);            # early withdrawal
    v = zeros(T, numSims);            # value
    a = zeros(T + 1,numSims);         # this is the path at the start of each period, so we include the 'start' of death
    b = zeros(T + 1,numSims);         # this is the path at the start of each period, so we include the 'start' of death

    # Other arrays that will be used below
    e     = zeros(T, numSims);        # the innovations to log income
    logZ1 = zeros(1, numSims);        # draws for initial income
    lZ    = zeros(Tretire, numSims)       # log income
    Z     = zeros(Tretire, numSims)       # stochastic part of income
    y    = zeros(T, numSims);        # income
    ypathIndex = zeros(T, numSims);   # holds the index (location) in the vector

    ## ------------------------------------------------------------------------
    # Obtain time series of incomes for our simulated individuals
    #-------------------------------------------------------------------------##
    # Draw random draws for starting income and for innovations
    seed1 = 1223424; # For the innovations
    seed2 = 234636;  # For initial income
    sig_inc = compute_sig_inc(params)
    e       = getNormalDraws( 0.0, sigma, T, numSims, seed1);    # normally distributed random draws for the innovation # NOTE: we actually just need Tretire
    logZ1   = getNormalDraws( mu, sigma, 1, numSims, seed2); # a random draw for the initial income
    # Later might replace sig_inc with 

    # Get all the incomes, recursively
    for s = 1:1: numSims                           # loop through individuals
        lZ[1, s] = truncate_custom(logZ1[1, s], mu - normBnd*sig_inc, mu + normBnd*sig_inc )
        Z[1, s]  = exp(lZ[1, s])

        for t = 1:1:Tretire-1                              # loop through time periods for a particular individual
            if (t != T)  # Get next year's income
                # ly[t+1, s] = (1 -rho) * mu + rho * ly[t, s] + e[t + 1, s];
                # ly[t+1, s] = truncate_custom(ly[t+1, s], -normBnd*sig_inc,normBnd*sig_inc );
                # y[t+1, s] = exp( ly[t+1, s] );

                lZ[t+1, s] = (1 -rho) * mu + rho * lZ[t, s] + e[t + 1, s]; # Note: the (1-rho) allows the AR(1) process to have a nonzero mean: http://www.maths.qmul.ac.uk/~bb/TimeSeries/TS_Chapter4_5.pdf
                lZ[t+1, s] = truncate_custom(lZ[t+1, s], mu - normBnd*sig_inc, mu + normBnd*sig_inc );  # NOTE: Cormac's code was missing a mu here. His version resulted in 0.0 starting income always
                Z[t+1, s] = exp( lZ[t+1, s] );
            end # if (t ~= T)

            if (t >= Tretire) 
                Z[t, s] = params["Yretire"]
            end # if (t >= Tretire)

            # Bound income if not extrapolating
            if !extrap_sim
                Z[t, s] = truncate_custom(Z[t, s], Ygrid[t,1], Ygrid[t,end])
            end
       end # t
    end # s

    # Now scale up based on deterministic component of income
    y[1:Tretire, :] = Z .* repeat(det_income, 1, numSims)

    for t = Tretire:T                              # loop through time periods for a particular individual
        y[t, :] = params["Yretire"] * ones(1, numSims)
    end

    ## ------------------------------------------------------------------------
    # Obtain consumption, asset and value profiles
    #-------------------------------------------------------------------------##

    # Initial values
    for s = 1:1: numSims
        a[1, s] = startA;
        b[1, s] = 0.0; 
    end

    for t = 1:1:T

        # (1) Define interpolation / extrapolation functions
        if (t < Tretire)                       

            # Define interpolation functions
            tA1 = policyA1[t, :, :, :];        # the relevant part of the liquid policy function
            tB1 = policyB1[t, :, :, :];        # the relevant part of the illiquid policy function
            tV  = EV[t, :, :, :];              # the relevant part of the value function

            knots_x = (Agrid[t, :], Bgrid[t, :], Ygrid[t, :])
            
            if extrap_sim
                # interpolate policy function, instructing to linearly extrapolate outside of bounds (can happen for income)            
                itp_A = interpolate(knots_x, tA1, Gridded(Linear()))
                itp_B = interpolate(knots_x, tB1, Gridded(Linear()))
                itp_V = interpolate(knots_x, tV, Gridded(Linear()))

                policyA1_fcn = extrapolate(itp_A,Line())
                policyB1_fcn = extrapolate(itp_B,Line())
                EV_fcn = extrapolate(itp_V,Line())
            else 
                # interpolate only
                policyA1_fcn = interpolate(knots_x, tA1, Gridded(Linear()))
                policyB1_fcn = interpolate(knots_x, tB1, Gridded(Linear()))
                EV_fcn       = interpolate(knots_x, tV, Gridded(Linear()))
            end
        else 
            tA1 = policyA1[t, :, :, 1];  # the relevant part of the policy function
            tB1 = policyB1[t, :, :, 1];  # the relevant part of the policy function
            tV  = EV[t, :, :, 1];        # the relevant part of the value function

            knots_x = (Agrid[t, :],Bgrid[t,:],)

            if extrap_sim
                itp_A = interpolate(knots_x, tA1, Gridded(Linear()))
                itp_B = interpolate(knots_x, tB1, Gridded(Linear()))
                itp_V = interpolate(knots_x, tV, Gridded(Linear()))

                policyA1_fcn = extrapolate(itp_A,Line())
                policyB1_fcn = extrapolate(itp_B,Line())
                EV_fcn       = extrapolate(itp_V,Line())
            else 
                policyA1_fcn = interpolate(knots_x, tA1, Gridded(Linear()))
                policyB1_fcn = interpolate(knots_x, tB1, Gridded(Linear()))
                EV_fcn       = interpolate(knots_x, tV, Gridded(Linear()))
            end
        end 

        # (2) Loop over all individuals using the above functions
        for s = 1:1: numSims
            if (t < Tretire) 
                # Working-Age Policy Functions - depend on income
                a[t+1, s] = policyA1_fcn(a[t, s], b[t,s], y[t, s])
                b[t+1, s] = policyB1_fcn(a[t, s], b[t,s], y[t, s])
                v[t, s] = EV_fcn(a[t, s], b[t,s], y[t, s])

            else
                # Retired Policy Functions - don't depend on income
                a[t+1, s] = policyA1_fcn(a[t, s], b[t,s])
                b[t+1, s] = policyB1_fcn(a[t, s], b[t,s])
                v[t, s]   = EV_fcn(a[t, s], b[t,s])
            end # if (t < Tretire)

            # Check whether next period's asset is below the lowest
            # permissable
            if ( a[t+1, s] < Agrid[t+1, 1] )
                println("Found a = $(a[t+1, s]) < 0 at time $t and sim $s given $(y[t,s])")
                a[t+1, s] = checkSimExtrap( Ygrid[t,:], Agrid[t+1, 1], y[t, s] );
            end

            # Impose constraint that illiquid assets must always be >= 0
            if ( b[t+1, s] < Bgrid[t+1, 1] )
                # println("Found b = $(b[t+1, s]) < 0 at time $t and sim $s given $(y[t,s])")
                b[t+1, s] = checkSimExtrap( Ygrid[t,:], Bgrid[t+1, 1], y[t, s] );
            end

            # Impose illiquid asset contribution limit
            # if ( b[t+1, s] > b[t, s] + max_contrib) 
            #     # println("b increase by more than max contrib")
            #     # println([ t, b[t+1, s] , b[t, s] , max_contrib ])
            #     b[t+1, s] = b[t, s] + max_contrib
            # end

            # Get consumption from today's assets, today's income and tomorrow's optimal assets
            # c[t, s] = a[t, s]  + y[t, s] - (a[t+1, s]/(1+r)) + b[t, s] - (b[t+1, s]/(1+params["r_b"])) - transaction_costs(params, t, Tretire, b[t+1, s], b[t, s]) 
            c[t, s] = (1.0+r)*a[t, s]  + y[t, s] - a[t+1, s] + (1+params["r_b"])*b[t, s] - b[t+1, s] - transaction_costs(params, t, Tretire, b[t+1, s], b[t, s]) 
            # NOTE: is it better to get the returns at the start or end of the period? Not sure... 
            
            # Record if early withdrawal
            # Note: could alternaively record if transaction_cost is positive
            B0 = b[t, s]
            B1 = b[t+1, s]
            B1_default = B0 * params["R_b"]
            if isapprox(B1_default, B1) | (t >= Tretire) 
                
            elseif (B1 < B1_default)
                ew[t,s] = 1.0
                cost = transaction_costs(params, t, Tretire, b[t+1, s], b[t, s])
                # println(cost)
                if cost <= 0.0
                    if (params["adj_cost_fixed"] > 0.0) | (params["adj_cost_prop"] > 0.0)
                        error("why no cost for early withdrawal?")
                    end
                end
            end

        end # s
    end #t

    return c, a, b, v, y, ew
end

function getNormalDraws( mu, sigma,  dim1, dim2, seed)
    # This function returns a dim1 * dim2 two array of pseudo random draws from
    # a normal distribution with mean mu and standard deviation sigma. A seed
    # is also inputted as well as lower and upper truncation points

    ## ------------------------------------------------------------------------
    #Set the seed
    seed!(seed)

    ## ------------------------------------------------------------------------
    # Draw standard normal draws, and transformthem so they come from a
    # distribution with our desired mean and stdev
    StdRandNormal = randn(dim1, dim2);
    normalDraws = mu .* ones(dim1, dim2)  + sigma .* StdRandNormal;
    return normalDraws
end

function checkSimExtrap( Ygrid, lba1, y )

    # This is a function to check, in selecting next period's asset we haven't
    # selected a value that is less than the borrowing constraint. This could
    # occur for one of two reasons. First it could be that we are extrapolating
    # (i.e. income in the period is larger than the largest income in the grid.
    # If this is the case we set next period's assets to the lowest permissable
    # level. Otherwise it is likely that there is an error - and in this case
    # we cause the programme to stop

    if (y > Ygrid[end]) || (y < Ygrid[1])
        a1 = lba1;
    else
        println( lba1, y )
        error("Next periods asset is below minimum permissable assets. And we are not extrapolating. Likely there is a bug")
    end
    return a1
end

function truncate_custom(y, negtrunc, postrunc)
    # Truncate input if its too big or too small
    if (y < negtrunc) # If y is less than the value negtrunc
       truncy = negtrunc;
    elseif (y > postrunc)
       truncy= postrunc; # If y is greater than the value postrunc
    else
       truncy = y;
    end
    return truncy
end
