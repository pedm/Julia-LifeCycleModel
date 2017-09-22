
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
        grid= linspace(minongrid, span, GridPoints)
    elseif method == "5logsteps"
        loggrid = linspace(log(1+log(1+log(1+log(1+log(1))))), log(1+log(1+log(1+log(1+log(1+span))))), GridPoints)
        grid = exp.(exp.(exp.(exp.(exp.(loggrid)-1)-1)-1)-1)-1
    end
end
