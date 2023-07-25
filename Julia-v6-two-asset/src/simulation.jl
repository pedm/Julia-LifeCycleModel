
function simNoUncer(params, Agrid, Ygrid, policyA1,V)

    ## ------------------------------------------------------------------------
    # Define parameters
    r       = params["r"]
    startA  = params["startA"]

    # Initialise arrays that will hold the paths of income consumption, value
    # and assets

    # Arguments for output
    c = zeros(T, 1)            # consumption
    v = zeros(T, 1)            # value
    a = zeros(T + 1,1)         # this is the path at the start of each period, so we include the 'start' of death
    y = zeros(T, 1)            # income

    ## ------------------------------------------------------------------------
    # Obtain paths using the initial condition and the policy and value
    # functions
    #-------------------------------------------------------------------------#
    a[1, 1] = startA
    for t = 1:1:T                     # loop through time periods for a particular individual
        
        # interpolate value function
        knots_x = (Agrid[t, :],)
        itp = interpolate(knots_x, V[t, :], Gridded(Linear()))
        v[t  , 1]   = itp[ a[t, 1] ]

        # interpolate the policy function
        knots_x = (Agrid[t, :],)
        itp_assets  = interpolate(knots_x, policyA1[t, :], Gridded(Linear()))
        a[t+1, 1]   = itp_assets[ a[t, 1] ]

        # fill in income and consumption
        y[t,   1]   = Ygrid[t];
        c[t  , 1]   = a[t, 1] + y[t, 1] - (a[t+1, 1]/(1+r))
    end   #t
    return c, a, v, y
end

function simWithUncer(params, Agrid, Bgrid, Ygrid, policyA1, EV)
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

    ## ------------------------------------------------------------------------
    # Initialise arrays that will hold the paths of income consumption, value
    # and assets

    # Arguments for output
    y = zeros(T, numSims);            # income
    c = zeros(T, numSims);            # consumption
    v = zeros(T, numSims);            # value
    a = zeros(T + 1,numSims);         # this is the path at the start of each period, so we include the 'start' of death
    b = zeros(T + 1,numSims);         # this is the path at the start of each period, so we include the 'start' of death

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
    sig_inc = sigma/((1-rho^2)^0.5); # long-run variance for AR(1), which is assumed variance of unconditional distribution at each age
    e = getNormalDraws( 0, sigma,  T, numSims, seed1);  # normally distributed random draws for the innovation
    logy1 =  getNormalDraws( mu, sig_inc,  1, numSims, seed2); # a random draw for the initial income

    # Get all the incomes, recursively
    for s = 1:1: numSims                           # loop through individuals
        ly[1, s] = truncate_custom(logy1[1, s], -normBnd*sig_inc,normBnd*sig_inc );
        y[1, s] = exp(ly[1, s]);
        for t = 1:1:T                              # loop through time periods for a particular individual
            if (t != T)  # Get next year's income
                ly[t+1, s] = (1 -rho) * mu + rho * ly[t, s] + e[t + 1, s];
                ly[t+1, s] = truncate_custom(ly[t+1, s], -normBnd*sig_inc,normBnd*sig_inc );
                y[t+1, s] = exp( ly[t+1, s] );
            end # if (t ~= T)

            if (t >= Tretire) # set income to zero if the individual has retired
                    y[t, s] = 0
            end # if (t >= Tretire)
       end # t
    end # s

## ------------------------------------------------------------------------
# Obtain consumption, asset and value profiles
#-------------------------------------------------------------------------##
     for s = 1:1: numSims
        a[1, s] = startA;
        b[1, s] = 1.0; # for some odd reason, it doesnt work if i set to zero.... b/c Bgrid doesn't go all the way down to zero, i think.

         for t = 1:1:T                              # loop through time periods for a particular individual

            if (t < Tretire)                       # first for the before retirement periods
                tA1 = policyA1[t, :, :, :];  # the relevant part of the policy function
                tV  = EV[t, :, :, :];         # the relevant part of the value function

                knots_x = (Agrid[t, :], Bgrid[t, :], Ygrid[t, :])
                
                # interpolate policy function, instructing to linearly extrapolate outside of bounds (can happen for income)
                itp = interpolate(knots_x, tA1, Gridded(Linear()))
                etp = extrapolate(itp,Line())
                a[t+1, s] = etp(a[t, s], b[t,s], y[t, s])

                # interpolate value function
                itp_V = interpolate(knots_x, tV, Gridded(Linear()))
                etp_V = extrapolate(itp_V,Line())
                v[t, s] = etp_V(a[t, s], b[t,s], y[t, s])

            else                          # next for the post retirement periods
                tA1 = policyA1[t, :, :, 1];  # the relevant part of the policy function
                tV = EV[t, :, :, 1];         # the relevant part of the value function

                knots_x = (Agrid[t, :],Bgrid[t,:],)

                # interpolate policy function
                itp = interpolate(knots_x, tA1, Gridded(Linear()))
                a[t+1, s] = itp(a[t, s], b[t,s])
                
                # interpolate value function
                itp_V = interpolate(knots_x, tV, Gridded(Linear()))
                v[t, s] = itp_V(a[t, s], b[t,s])

            end # if (t < Tretire)

            # Check whether next period's asset is below the lowest
            # permissable
            if ( a[t+1, s] < Agrid[t+1, 1] )
               a[t+1, s] = checkSimExtrap( Ygrid, Agrid[t+1, 1], y[t, s] );
            end

            # Get consumption from today's assets, today's income and
            # tomorrow's optimal assets
            c[t, s] = a[t, s]  + y[t, s] - (a[t+1, s]/(1+r));
            b[t+1, s] = b[t, s] # TODO: make this later...
        end   #t
    end # s

    return c, a, b, v, y
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
