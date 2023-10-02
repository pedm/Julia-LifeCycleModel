# TODO: will there be a speed improvement if I specify that params is Dict{Any,Any} ?
@everywhere function utility(params::Dict{String,Float64}, cons::Float64)
    # Note: seems we save a lot of time by specifying type for params

    if params["gamma"] == 1
         utils = log(cons)
    else
        utils = ((cons).^(params["gamma_mod"])  )./(params["gamma_mod"])
    end
end

@everywhere function objectivefunc(params::Dict{String,Float64}, itp, A1::Float64, A0::Float64, Y::Float64)

    #-------------------------------------------------------------------------------#
    # This function returns the following quantity:
    # - (u(c) +  b V( A1))
    # where c is calculated from today's assets and tomorrow's assets

    #Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    #total value (u(c) + b * VA1
    # beta = params["beta"]
    # r = params["r"]

    cons = A0 + Y - (A1)/(1 + params["r"])
    value = utility(params, cons) + params["beta"] * itp[A1]

    ## ------------------------------------------------------------------------
    #The optimisation routine that we will use searches for the minimum of the
    #function. We want the maximum. So we multiply out function here by -1 so
    #that the optimiser will fill the minimum of the negative of our function,
    #i.e. the maximum of our functino

    value = - value
    return value
end

@everywhere function transaction_costs(B1, B0)

    # Most realistic: Only have a fee if withdrawing:
    # TODO: maybe later: want to have some subsidy to contributions
    if isapprox(B0, B1) 
        return 0.0
    elseif B1 < B0
        return params["adj_cost_fixed"] + params["adj_cost_prop"] * abs(B1 - B0) # TODO: might not need that abs()
    else
        return 0.0
    end
end
