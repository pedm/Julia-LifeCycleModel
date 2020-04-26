
# TODO: add params
function simNoUncer(Agrid, Ygrid, policyA1,V,startA)
    # Initialise arrays that will hold the paths of income consumption, value
    # and assets

    # Arguments for output
    c = zeros(T, 1)            # consumption
    v = zeros(T, 1)            # value
    a = zeros(T + 1,1)         # this is the path at the start of each period, so we include the 'start' of death
    y = zeros(T, 1)              # income

    ## ------------------------------------------------------------------------
    # Obtain paths using the initial condition and the policy and value
    # functions
    #-------------------------------------------------------------------------#
    a[1, 1] = startA
    for t = 1:1:T                     # loop through time periods for a particular individual
        # Original matlab code
        # v[t  , 1]   = interp1(Agrid[t, :],V[t, :],a[t, 1],interpMethod, 'extrap')
        knots_x = (Agrid[t, :],)
        itp = interpolate(knots_x, V[t, :], Gridded(Linear()))
        v[t  , 1]   = itp( a[t, 1] )

        # a[t+1, 1]   = interp1(Agrid[t, :], policyA1[t, :] ,a[t, 1],interpMethod, 'extrap')
        knots_x = (Agrid[t, :],)
        itp_assets  = interpolate(knots_x, policyA1[t, :], Gridded(Linear()))
        a[t+1, 1]   = itp_assets[ a[t, 1] ]
        y[t,   1]   = Ygrid[t];
        c[t  , 1]   = a[t, 1] + y[t, 1] - (a[t+1, 1]/(1+r))
    end   #t
    return c, a, v, y
end

function simWithUncer(params, Agrid, Ygrid, policyA1, EV)
    # (policyA1,EV,startA)

    # This function takes the policy functions and value functions, along with
    # starting assets and returns simulated paths of income, consumption,
    # assets and value

    ## ------------------------------------------------------------------------
    # Define parameters

    rho     = params["rho"]
    mu      = params["mu"]
    sigma   = params["sigma"]
    r       = params["r"]
    startA  = params["startA"]

    # Setup nodes and weights for transitory income shocks
    μtransshocks = -0.5 * params["sigma_trans"]
    Ytrans_grid, Ytrans_weights = gausshermite_normal_distribution(numPointsYTrans, μtransshocks, params["sigma_trans"] )
    Ytrans_grid = exp.(Ytrans_grid)

    ## ------------------------------------------------------------------------
    # Initialise arrays that will hold the paths of income consumption, value
    # and assets

    # Arguments for output
    y = zeros(T, numSims);            # income
    y_perm = zeros(T, numSims);            # persistent component of income
    c = zeros(T, numSims);            # consumption
    v = zeros(T, numSims);            # value
    a = zeros(T + 1,numSims);         # this is the path at the start of each period, so we include the 'start' of death

    # Other arrays that will be used below
    e     = zeros(T, numSims);        # the innovations to log income
    logy1 = zeros(1, numSims);        # draws for initial income
    ly    = zeros(T, numSims);        # log income
    ypathIndex = zeros(T, numSims);   # holds the index (location) in the vector

    ## ------------------------------------------------------------------------
    # Obtain time series of incomes for our simulated individuals
    #-------------------------------------------------------------------------##
    # Draw random draws for starting income and for innovations
    seed1 = 1223424; # For the innovations
    seed2 = 234636;  # For initial income
    seed3 = 2354654;  # For transitory income shocks

    sig_inc        = sigma/ ((1-rho^2)^0.5)
    e              = getNormalDraws( 0, sigma, T, numSims, seed1); # normally distributed random draws for the innovation
    logy1          = getNormalDraws( mu, sig_inc, 1, numSims, seed2); # a random draw for the initial income

    ly_trans       = getNormalDraws( μtransshocks, params["sigma_trans"], T, numSims, seed2) # a random draw for the transitory income shocks
    ly_trans_trunc = map(x -> truncate_custom(x, μtransshocks-normBnd*params["sigma_trans"], μtransshocks+normBnd*params["sigma_trans"] ), ly_trans)
    y_trans        = exp.(ly_trans_trunc)

    # Get all the incomes, recursively
    for s = 1:1: numSims                           # loop through individuals
        ly[1, s] = truncate_custom(logy1[1, s], -normBnd*sig_inc,normBnd*sig_inc );
        y_perm[1, s] = exp(ly[1, s]);
        for t = 1:1:T                              # loop through time periods for a particular individual
            if (t != T)  # Get next year's income
                ly[t+1, s]     = (1 -rho) * mu + rho * ly[t, s] + e[t + 1, s];
                ly[t+1, s]     = truncate_custom(ly[t+1, s], -normBnd*sig_inc,normBnd*sig_inc );
                y_perm[t+1, s] = exp( ly[t+1, s] )
            end # if (t ~= T)

            if (t >= Tretire) # set income to zero if the individual has retired
                    y_perm[t, s] = 0.0
                    y_trans[t, s] = 0.0
            end # if (t >= Tretire)
       end # t
    end # s

    y = y_perm .* y_trans

## ------------------------------------------------------------------------
# Obtain consumption, asset and value profiles
#-------------------------------------------------------------------------##
     for s = 1:1: numSims
        a[1, s] = startA;
         for t = 1:1:T                              # loop through time periods for a particular individual
            if (t < Tretire)                       # first for the before retirement periods
                # clear tA1 tV;                      #necessary as the dimensions of these change as we wor through this file
                tA1 = policyA1[t, :, :, :];  # the relevant part of the policy function
                tV  = EV[t, :, :];         # the relevant part of the value function
                # a[t+1, s] =  interp2D(Agrid[t,:]', Ygrid[t, :]', tA1, a[t, s], y[t, s]);
                # v[t  , s] =  interp2D(Agrid[t,:]', Ygrid[t, :]', tV , a[t, s], y[t, s]);

                knots_x = (Agrid[t, :], Ygrid[t, :], Ytrans_grid)
                knots_EV = (Agrid[t, :], Ygrid[t, :])

                itp_basic = interpolate(knots_x, tA1, Gridded(Linear()))
                itp       = extrapolate(itp_basic, Interpolations.Line() ) # Extrapolation necessary because Ygrid points are in the middle of the two extremes

                try
                    a[t+1, s] = itp( a[t, s], y_perm[t, s], y_trans[t, s] )
                catch
                    println("Some issue with interpolation: print diagonistics")
                    println("t = $t              ")
                    println("x = $(knots_x      )")
                    println("a = $(a[t, s]      )")
                    println("yperm = $(y_perm[t, s] )")
                    println("y trans = $(y_trans[t, s])")

                    # NOTE: April 2020: seems like the issue is that yperm < the smallest grid point. WHY !?!?!
                    # TODO: figure out what causes that

                    a[t+1, s] = itp( a[t, s], y_perm[t, s], y_trans[t, s] )
                end

                itp_V_basic = interpolate(knots_EV, tV, Gridded(Linear()))
                itp_V       = extrapolate(itp_V_basic, Interpolations.Flat() )

                v[t, s] = itp_V( a[t, s], y_perm[t, s] )

            else                          # next for the post retirement periods
                # clear tA1 tV;
                tA1 = policyA1[t, :, 1, 1];  # the relevant part of the policy function
                tV = EV[t, :, 1];         # the relevant part of the value function
                # a[t+1, s] = interp1(Agrid[t,:]', tA1, a[t, s], 'linear', 'extrap');
                # v[t,   s] = interp1(Agrid[t,:]', tV , a[t, s], 'linear', 'extrap');

                knots_x = (Agrid[t, :],)
                itp = interpolate(knots_x, tA1, Gridded(Linear()))
                a[t+1, s] = itp( a[t, s] )

                itp_V = interpolate(knots_x, tV, Gridded(Linear()))
                v[t, s] = itp_V( a[t, s] )

            end # if (t < Tretire)

            # Check whether next period's asset is below the lowest
            # permissable
            if ( a[t+1, s] < Agrid[t+1, 1] )
               a[t+1, s] = checkSimExtrap( Ygrid, Agrid[t+1, 1], y[t, s] );
            end

            # Get consumption from today's assets, today's income and
            # tomorrow's optimal assets
            c[t, s] = a[t, s]  + y[t, s] - (a[t+1, s]/(1+r));
        end   #t
    end # s

    return c, a, v, y
end

function getNormalDraws( mu, std_dev,  dim1, dim2, seed)
    # This function returns a dim1 * dim2 two array of pseudo random draws from
    # a normal distribution with mean mu and standard deviation sigma. A seed
    # is also inputted as well as lower and upper truncation points

    ## ------------------------------------------------------------------------
    #Set the seed
    # srand(seed)
    Random.seed!(seed)

    ## ------------------------------------------------------------------------
    # Draw standard normal draws, and transformthem so they come from a
    # distribution with our desired mean and stdev
    StdRandNormal = randn(dim1, dim2);
    normalDraws = mu .* ones(dim1, dim2)  + std_dev .* StdRandNormal;
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

    if (y > Ygrid[numPointsY]) || (y < Ygrid[1])
        a1 = lba1;
    else
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
