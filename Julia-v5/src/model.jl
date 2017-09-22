function utility(cons::Float64)
    # if gamma == 1
    #     utils = log(cons)
    # else
        utils = ((cons)^(gamma_mod)  )/(gamma_mod)
    # end
end


function objectivefunc(itp, A1, A0, Y)

    #-------------------------------------------------------------------------------#
    # This function returns the following quantity:
    # - (u(c) +  b V( A1))
    # where c is calculated from today's assets and tomorrow's assets

    #Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    #total value (u(c) + b * VA1
    # beta = params["beta"]
    # r = params["r"]

    cons = A0 + Y - (A1)/(1+r)
    value = utility(cons) + beta * itp[A1]

    ## ------------------------------------------------------------------------
    #The optimisation routine that we will use searches for the minimum of the
    #function. We want the maximum. So we multiply out function here by -1 so
    #that the optimiser will fill the minimum of the negative of our function,
    #i.e. the maximum of our functino

    value = - value
    return value
end
