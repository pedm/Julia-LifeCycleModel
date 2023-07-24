

function simNoUncer(T, r, Agrid, policyA1,V,startingA)
    # Initialise arrays that will hold the paths of income consumption, value
    # and assets

    # Arguments for output
    c = zeros(T, 1)            # consumption
    v = zeros(T, 1)            # value
    a = zeros(T + 1,1)         # this is the path at the start of each period, so we include the 'start' of death

    ## ------------------------------------------------------------------------
    # Obtain paths using the initial condition and the policy and value
    # functions
    #-------------------------------------------------------------------------#
    a[1, 1] = startingA
    for t = 1:1:T                     # loop through time periods for a particular individual
        # Original matlab code
        # v[t  , 1]   = interp1(Agrid[t, :],V[t, :],a[t, 1],interpMethod, 'extrap')
        knots_x = (Agrid[t, :],)
        itp = interpolate(knots_x, V[t, :], Gridded(Linear()))
        v[t  , 1]   = itp(a[t, 1])

        # a[t+1, 1]   = interp1(Agrid[t, :], policyA1[t, :] ,a[t, 1],interpMethod, 'extrap')
        knots_x = (Agrid[t, :],)
        itp_assets  = interpolate(knots_x, policyA1[t, :], Gridded(Linear()))
        a[t+1, 1]   = itp_assets(a[t, 1])

        c[t  , 1] = a[t, 1]  - (a[t+1, 1]/(1+r))
    end   #t
    return c, a, v
end
