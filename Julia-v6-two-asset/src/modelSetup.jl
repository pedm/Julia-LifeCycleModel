
function getMinAndMaxAss(model, minInc, maxInc)

    params     = model["params"]
    # ints       = model["ints"]
    # T          = ints["T"]

    ## ------------------------------------------------------------------------
    # Initialise the output matrices
       BC   = zeros(T+1,1)
       maxA = zeros(T+1,1)
       maxB = zeros(T+1,1)

    ## ------------------------------------------------------------------------
    # Iteratively, calculate the borrowing constraints, and maximum on asset
    # grid

    # Borrowing constraints
    BC[T + 1] = 0
    for ixt = T:-1:1
        BC[ixt] = BC[ixt+1]/(1 + params["r"]) - minInc[ixt, 1] + params["minCons"]
    end

    # if borrowing is not allowed, replace negative points in the borrowing
    # constraint with zero
    if (borrowingAllowed == 0 )
        BC[BC.<0] .= 0
    end

    # HARD CODE ZEROS for BC = minAss
    # BC = zeros(size(BC))

    # Maximum illiquid assets
    # maxB[1] = 0.1 # arbitrary. just have this to prevent the first period grid being all zeros
    # for ixt = 2:1:T+1
    #     maxB[ixt] = (maxB[ixt - 1] + params["max_contrib"] ) * (1 + params["r_b"])      
    # end

    # if params["r"] > params["r_b"]
    #     error("current grid setup won't support these params")
    # end

    # Maximum Liquid Assets
    maxA[1] = params["startA"]
    # If starting assets are 0 we will have maxA<=0 in period 1 and BC>0.
    # Therefore BC>maxA (minA<maxA). Ensure that this is not the case by
    # replacing any element of maxA that is less than the corresponding element
    # of BC with that corresponding element plus 1
    maxA[1] = params["startA"] 
    if maxA[1] <= BC[1]
        maxA[1] = BC[1] + 1
    end

    max_r = maximum([params["r"], params["r_b"]])

    for ixt = 2:1:T+1
        # Original version: one asset
        # maxA[ixt] = (maxA[ixt - 1] + maxInc[ixt-1] ) * (1 + params["r"])      

        # Updated version: two asset, one with higher return (assume that the agent return maximizes, so makes maximum contribution to liquid... 
        # then later add on maxB to capture the possibility that they withdraw all that period)
        maxA[ixt] = (maxA[ixt - 1] + maxInc[ixt-1]) * (1.0 + max_r)      

        if maxA[ixt] <= BC[ixt]
          maxA[ixt] = BC[ixt] + 1
        end
    end

    # Now allow for the possiblity that you extract everything from their retirment account to liquid account (again all under assumption that you get higher returns on retirement account)
    # maxA = maxA .+ maxB
    maxB = maxA

    return BC, maxA, maxB
end

function getGrid(minongrid, maxongrid, GridPoints, method)
    span = maxongrid - minongrid
    if method == "equalsteps"
        grid= collect(range(minongrid, maxongrid, GridPoints))
    elseif method == "logsteps"
        loggrid = collect(range(log(1) , log(1+span), GridPoints))
        grid = minongrid.+ exp.(loggrid) .-1
    elseif method == "5logsteps"
        loggrid = collect(range(log(1+log(1+log(1+log(1+log(1))))), log(1+log(1+log(1+log(1+log(1+span))))), GridPoints))
        grid = minongrid.+ exp.(exp.(exp.(exp.(exp.(loggrid) .-1) .-1) .-1) .-1) .-1
    end
end

function tauchen(params)

    # First get the standard deviation of income (from sigma and rho) (note: this is standard deviation of income in the steady state)
    # We confirmed in Adda & Cooper that this is indeed std and not variance
    sig_inc = compute_sig_inc( params )
    # NOTE: sig_inc will be Inf if rho = 1.
    # sig_inc is used to define Ygrid

    # Split the entire normal distribution into numPointsY sections that
    # are equiprobable. The output lNormDev gives the (numPointsY + 1)
    # points that bound the sections, the output ly gives the
    # (numPointsY) expected value in each section
    lNormDev, lZ = getNormDev(params["mu"], sig_inc, normBnd, numPointsY )

    #---------------------#
    #Get transition matrix Q(i, j). The prob of income j in t+1
    #conditional on income i in t
    #---------------------#
    Q = zeros(numPointsY, numPointsY)              #initialise the transition matrix
    for i = 1:1:numPointsY
        for j = 1:1:numPointsY
            hiDraw = lNormDev[j+1] - (1-params["rho"])*params["mu"] - params["rho"] * lZ[i]; #highest innovation that will give us income j tomorrow
            loDraw = lNormDev[j]   - (1-params["rho"])*params["mu"] - params["rho"] * lZ[i]; #lowest  innovation that will give us income j tomorrow
            Q[i,j] = stdnormcdf_manual(hiDraw/params["sigma"]) - stdnormcdf_manual(loDraw/params["sigma"]);
        end #j

        #Each of the rows of Q should add up to 1. But
        #due to the truncation of the normal distribution they will
        #not. So we divide through by the sum of elements in the row
        Q[i, :] = Q[i, :] ./ sum(Q[i, :])
    end #i

    return lZ, Q 
end

function getIncomeGrid(model)

    # A function that returns:
    # 1. an income grid
    # 2. A Markovian transition matrix (Q) over income realisations
    # 3. A vector of minimum incomes in each year
    # 4. A vector of maximum incomes in each year

    params     = model["params"]

    #----------------------------------------#
    # Income draws are log normally distributed
    #----------------------------------------#
    sig_inc = compute_sig_inc(params)

    # lZ, Q = tauchen(params)
    lZ, Q = rouwenhorst(numPointsY, params["rho"], params["sigma"])

    Z = exp.(lZ);                   # Get y from log y
    minZ = exp(-normBnd * sig_inc); # Get the minimum income in each year # NOTE: are we missing mu here? Cormac bug ??
    maxZ = exp(normBnd * sig_inc);  # Get the maximum income in each year

    if (Z[1] < 1e-4) || (Z[ numPointsY ] > 1e5)
        warn("Combination of sigma and rho give a very high income variance. Numerical instability possible", once = true)
    end

    #----------------------------------------#
    #Now get a matrix, T * numPointsY that holds the grid for each income in
    #each year. Do likewise with minimum and maximum income
    #----------------------------------------#
    Zgrid = repeat(Z',      Tretire)
    minZ  = repeat([minZ'], Tretire)
    maxZ  = repeat([maxZ'], Tretire)


    #----------------------------------------#
    # Add deterministic age trend
    #----------------------------------------#

    # Estimated on PSID using ages 22 - 65 
    # from PSID i get: logincome = 8.247418 + .1123966 * Age + -.0114392  * (Age^2)/10
    if Tretire == 44
        t         = collect(22:1:22+Tretire - 1 )  # this works when Tretire = 44
        logincome = params["inc_reg_constant"]*ones(size(t)) + params["inc_reg_age"] * t + (params["inc_reg_age2"] * t.^2) + (params["inc_reg_age3"] * t.^3) # + (params["inc_reg_age4"] * t.^4)/1000
    else
        # this will be the same as the above case when Tretire = 44
        t = collect(LinRange(22, 20+44-1 , Tretire )) # get evenly spaced t
        logincome = params["inc_reg_constant"]*ones(size(t)) + params["inc_reg_age"] * t + (params["inc_reg_age2"] * t.^2) + (params["inc_reg_age3"] * t.^3) # + (params["inc_reg_age4"] * t.^4)/1000
    end

    income = exp.(logincome)

    # Store the value of initial income -- later we will scale up using this value
    income_scale = copy(income[1])
    println("income_scale = $income_scale")

    # Why do we do this? Answer: when we create 5logstep asset grid, very different grid if income is 25k per year vs 1-2 per year
    det_income = income ./ income_scale # standardized_income

    #----------------------------------------#
    #Now get a matrix, T * numPointsY that holds the grid for each income in
    #each year. Do likewise with minimum and maximum income
    #----------------------------------------#
    Ygrid = zeros(T, numPointsY)
    minInc = zeros(T)
    maxInc = zeros(T)

    Ygrid[1:Tretire, :]  = Zgrid .* repeat( det_income, 1, numPointsY)            # construct Ygrid pre ret
    minInc[1:Tretire]  = minZ .* det_income
    maxInc[1:Tretire]  = maxZ .* det_income
    # ret_income = params["ret_fraction"] * det_income[end] * ones(T - Tretire)

    #----------------------------------------#
    #Replace these arrays with zeros for all years after retirement
    #----------------------------------------#
    Ygrid[Tretire:T, :] .= params["Yretire"]
    minInc[Tretire:T] .= params["Yretire"]
    maxInc[Tretire:T] .= params["Yretire"]

    return Ygrid, Q, minInc, maxInc, det_income

end

function getAssetGrid(model)

    params     = model["params"]
    # ints       = model["ints"]
    # T          = ints["T"]
    # numPointsA = ints["numPointsA"]
    # numPointsB = ints["numPointsB"]
    
    MinAss, MaxAss, MaxB = getMinAndMaxAss(model, minInc, maxInc)
    Agrid = zeros(T+1, numPointsA)
    Bgrid = zeros(T+1, numPointsB)
    for ixt = 1:1:T+1
        Agrid[ixt, :] = getGrid(MinAss[ixt], MaxAss[ixt], numPointsA, gridMethod)
        Bgrid[ixt, :] = getGrid(MinAss[ixt], MaxB[ixt], numPointsB, gridMethod)
    end

    return Agrid, Bgrid
end

function compute_sig_inc(params::Dict{String,Float64})
    # We use sig_inc to define the grid points in Ygrid and to define the places at which we truncate the income distribution (aka we truncate income at normBnd*sig_inc)

    rho = params["rho"]
    sig_inc = params["sigma"] / (( 1.0 - rho ^2.0 )^0.5)
    return sig_inc
end

function getNormDev(mu, sigma_inc, trunc, N )

    # This function returns two vectors:
    # a) Z - a vector of (N+1) normal deviates (numbers) that divide a normal
    # distribution with standard deviation sigma into N segments,
    # each of which are equiprobable
    # b) EVbetweeenZ - the expected value of the random variables between those
    # two points


    ## -------------------------------------------------------------------------
    # Initialise the output and working arrays

    # Output
    Z = zeros(N + 1, 1)
    EVbetweeenZ = zeros(N, 1);

   ## -------------------------------------------------------------------------
   #Finding the points that divide the standard normal into N segments
   #The first and last of these should, if we are using an actual normal
   #distribution would be minus and plus infinity.  In reality we use a truncated
   #normal distribution and use -trunc and +trunc (where trunc is an argument)
   Z[1] = -trunc * sigma_inc;
   Z[N + 1] = trunc * sigma_inc;

    #Now recursively get the rest of the points
    for ixi = 2:1:N
        Z[ixi] = sigma_inc * stdnorminv_manual((ixi-1)/N) + mu;
    end

    stdZ = (Z - mu .* ones(N + 1, 1)) ./ sigma_inc;

    #Finding the expected value within each interval (see Adda & Cooper page
#58)
    PDF = stdnormpdf_manual(stdZ);

    for ixi = 1:1:N
        EVbetweeenZ[ixi] = N .* sigma_inc .* (PDF[ixi] - PDF[ixi + 1]) + mu
    end

    return Z, EVbetweeenZ
end

function stdnorminv_manual(p)
    # This routines gives the standard normal value associated with a particular p value,
    # with truncations at -3 and 3

    if (p<stdnormcdf_manual(-3))
        x = -3;
    elseif  (p>stdnormcdf_manual(3))
        x = 3;
    else
        # boundforzero = [-3, 3];
        # fzero  Single-variable nonlinear zero finding.
        # x = fzero(@(t) stdnormcdf_manual(t) - p, boundforzero, optimset('TolX',1e-12))
        f(x) = stdnormcdf_manual(x) - p

        # NOTE: for some reason, f(3) and f(3.0) gave different results back when using simps_custom(). thats why I moved to quadgk()
        # display(f(-3.0))
        # display(f(3.0))

        x = fzero(f, -3.0, 3.0)
    end
    return x
end

function stdnormcdf_manual(x)
    # This function gives the cdf of a standard normal
    approxforminusinf = -20;
    # quad   Numerically evaluate integral, adaptive Simpson quadrature.
    # cdf = quad(@(t) stdnormpdf_manual(t), approxforminusinf, x, 1e-12);
    # cdf = simps_custom(t -> stdnormpdf_manual(t), approxforminusinf, x, 1000)

    cdf, error = quadgk(stdnormpdf_manual, approxforminusinf, x)
    return cdf
end

function stdnormpdf_manual(x)
    # This function gives the pdf of a standard normal
    pdf = ((sqrt(2*pi))^-1) * (exp.( - ((x).^2)/(2)))
    return pdf
end

# Quadrature integration (replicate quad() from matlab in julia)
# Copied from http://blog.mmast.net/simpson-integration-julia
# Note: could also use https://lectures.quantecon.org/jl/julia_libraries.html
function simps_custom(f::Function, a::Number, b::Number, n::Number)
    # f = function over which you want to integrate between a and b
    # n = number of grid points

    # n # 2 == 0 || error("`n` must be even")
    h = (b-a)/n
    s = f(a) + f(b)
    s += 4sum(f(a + collect(1:2:n) * h))
    s += 2sum(f(a + collect(2:2:n-1) * h))
    return h/3 * s
end
