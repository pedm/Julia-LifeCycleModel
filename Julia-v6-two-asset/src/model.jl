# TODO: will there be a speed improvement if I specify that params is Dict{Any,Any} ?
@everywhere function utility(params::Dict{String,Float64}, cons::Float64)
    # Note: seems we save a lot of time by specifying type for params

    if params["gamma"] == 1
         utils = log(cons)
    else
        utils = ((cons).^(params["gamma_mod"])  )./(params["gamma_mod"])
    end
end

@everywhere function objectivefunc(params::Dict{String,Float64}, itp, A1::Float64, Atilde::Float64, B1::Float64)

    #-------------------------------------------------------------------------------#
    # This function returns the following quantity:
    # - (u(c) +  b V( A1))
    # where c is calculated from today's assets and tomorrow's assets


    cons  = Atilde - A1 
    value = utility(params, cons) + params["beta"] * itp[A1, B1]
    # value = utility(params, cons) + params["beta"] * itp[A1]

    ## ------------------------------------------------------------------------
    #The optimisation routine that we will use searches for the minimum of the
    #function. We want the maximum. So we multiply out function here by -1 so
    #that the optimiser will fill the minimum of the negative of our function,
    #i.e. the maximum of our function

    value = - value
    return value
end

@everywhere function transaction_costs(ixt, B1, B0)

    # Most realistic: Only have a fee if withdrawing:
    # TODO: maybe later: want to have some subsidy to contributions

    # TODO: should define B1_default
    B1_default = B0 * (1.0+params["r_b"])

    if isapprox(B1_default, B1) # | (ixt >= Tretire)
        return 0.0
    elseif B1 < B0 # B1_default
        early_withdrawal = abs(B1 - B1_default) # TODO: might not need that abs()
        return params["adj_cost_fixed"] + params["adj_cost_prop"] * early_withdrawal 
    else
        return 0.0
    end
end
