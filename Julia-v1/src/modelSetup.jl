
function getMinAndMaxAss()
    ## ------------------------------------------------------------------------
    # Initialise the output matrices
       BC   = zeros(T+1,1)
       maxA = zeros(T+1,1)
    # do I have to initialise the size NaN(T+1, 1) ?

    ## ------------------------------------------------------------------------
    # Iteratively, calculate the borrowing constraints, and maximum on asset
    # grid

    # Borrowing constraints
    BC[T + 1] = 0
    for ixt = T:-1:1
        BC[ixt] = BC[ixt+1]/(1+r) + minCons
    end

    # Maximum Assets
    maxA[1] = startA
    for ixt = 2:1:T+1
        maxA[ixt] = (maxA[ixt - 1] ) * (1+r)
    end
    return BC, maxA
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
