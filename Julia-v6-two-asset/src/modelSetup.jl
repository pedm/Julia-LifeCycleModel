
function getMinAndMaxAss(params, minInc, maxInc)

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


function getIncomeGrid(params)

    # A function that returns:
    # 1. an income grid
    # 2. A Markovian transition matrix (Q) over income realisations
    # 3. A vector of minimum incomes in each year
    # 4. A vector of maximum incomes in each year

    #----------------------------------------#
    # Scenario where there is no uncertainty
    #----------------------------------------#
    if isUncertainty == 0
        y = exp( params["mu"] )                 # income is set equal to the exp of the log mean
        minInc = y
        maxInc = y
        Q = [1]                    # The transition matrix Q is simply a constant 1
                                   # with prob 1 each period income is 1

       #----------------------------------------#
       #Now get a matrix, T * numPointsY that holds the grid for each income in
       #each year. Do likewise with minimum and maximum income
       #----------------------------------------#
       Ygrid = repeat([y'], T, 1)
       minInc = repeat([minInc'], T, 1)
       maxInc = repeat([maxInc'], T, 1)

    #----------------------------------------#
    # Scenario where there is uncertainty - income draws are log normally distributed
    #----------------------------------------#

elseif isUncertainty == 1

           # First get the standard deviation of income (from sigma and rho)
           sig_inc = params["sigma"] / ((1- params["rho"]^2)^0.5)

           # Split the entire normal distribution into numPointsY sections that
           # are equiprobable. The output lNormDev gives the (numPointsY + 1)
           # points that bound the sections, the output ly gives the
           # (numPointsY) expected value in each section
           lNormDev, ly = getNormDev(params["mu"], sig_inc, normBnd, numPointsY );

           #---------------------#
           #Get transition matrix Q(i, j). The prob of income j in t+1
           #conditional on income i in t
           #---------------------#
           Q = zeros(numPointsY, numPointsY)              #initialise the transition matrix
           for i = 1:1:numPointsY
               for j = 1:1:numPointsY
                   hiDraw = lNormDev[j+1] - (1-params["rho"])*params["mu"] - params["rho"] * ly[i]; #highest innovation that will give us income j tomorrow
                   loDraw = lNormDev[j]   - (1-params["rho"])*params["mu"] - params["rho"] * ly[i]; #lowest  innovation that will give us income j tomorrow
                   Q[i,j] = stdnormcdf_manual(hiDraw/params["sigma"]) - stdnormcdf_manual(loDraw/params["sigma"]);
               end #j

               #Each of the rows of Q should add up to 1. But
               #due to the truncation of the normal distribution they will
               #not. So we divide through by the sum of elements in the row
               Q[i, :] = Q[i, :] ./ sum(Q[i, :])
           end #i

           y = exp.(ly);                      # Get y from log y
           minInc = exp(-normBnd * sig_inc); #Get the minimum income in each year
           maxInc = exp(normBnd * sig_inc);  #Get the maximum income in each year

           if (y[1] < 1e-4) || (y[ numPointsY ] > 1e5)
               warning("Combination of sigma and rho give a very high income variance. Numerical instability possible")
           end

           #----------------------------------------#
           #Now get a matrix, T * numPointsY that holds the grid for each income in
           #each year. Do likewise with minimum and maximum income
           #----------------------------------------#
           Ygrid = repeat(y', T, 1)
           minInc = repeat([minInc], T, 1)
           maxInc = repeat([maxInc], T, 1)

    end  # if isUncertainty == 0



    #----------------------------------------#
    #Replace these arrays with zeros for all years after retirement
    #----------------------------------------#

    if Tretire == 0         # no work (retired at birth)
       Ygrid[:, :] .= params["Yretire"]
       minInc[:, :] .= params["Yretire"]
       maxInc[:, :] .= params["Yretire"]
    elseif (Tretire > 0) && (Tretire <=T)  #retire at some age
        Ygrid[Tretire:T, :] .= params["Yretire"]
        minInc[Tretire:T, :] .= params["Yretire"]
        maxInc[Tretire:T, :] .= params["Yretire"]
    end

    return Ygrid, Q, minInc, maxInc

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
