@everywhere function utility(params::Dict{String,Float64}, cons::Float64, cons_tempt::Float64)
    # utils = felicity(params, cons)

    if cons_tempt < cons
        error("Why is cons_tempt < cons?")
    end

    utils = (1.0 + params["lambda"]) * felicity(params, cons) - params["lambda"] * felicity(params, cons_tempt)
    return utils
end

@everywhere function felicity(params::Dict{String,Float64}, cons::Float64)

    if params["gamma"] == 1
         utils = log(cons)
    else
        utils = ((cons).^(params["gamma_mod"])  )./(params["gamma_mod"])
    end
end

@everywhere function objectivefunc(params::Dict{String,Float64}, itp, A1::Float64, Atilde::Float64, B1::Float64, cons_tempt::Float64)

    #-------------------------------------------------------------------------------#
    # This function returns the following quantity:
    # - (u(c) +  b V( A1))
    # where c is calculated from today's assets and tomorrow's assets


    cons  = Atilde - A1 
    value = utility(params, cons, cons_tempt) + params["beta"] * itp[A1, B1]
    # value = utility(params, cons) + params["beta"] * itp[A1, B1]
    # value = utility(params, cons) + params["beta"] * itp[A1]

    ## ------------------------------------------------------------------------
    #The optimisation routine that we will use searches for the minimum of the
    #function. We want the maximum. So we multiply out function here by -1 so
    #that the optimiser will fill the minimum of the negative of our function,
    #i.e. the maximum of our function

    value = - value
    return value
end

@everywhere function transaction_costs(params, ixt, Tretire, B1, B0)

    # Most realistic: Only have a fee if withdrawing:
    # TODO: maybe later: want to have some subsidy to contributions

    B1_default = B0 * params["R_b"]

    # Orig version: 
    # if isapprox(B1_default, B1)
    #     return 0.0
    # elseif B1 < B0 # B1_default
    #     early_withdrawal = abs(B1 - B1_default) # TODO: might not need that abs()
    #     return params["adj_cost_fixed"] + params["adj_cost_prop"] * early_withdrawal 
    # else
    #     return 0.0
    # end

    if isapprox(B1_default, B1) # | (ixt >= Tretire)
        return 0.0
    elseif B1 < B1_default
        early_withdrawal = abs(B1 - B1_default) # TODO: might not need that abs()

        if (ixt >= Tretire)
            return params["ret_adj_cost_fixed"] + params["ret_adj_cost_prop"] * early_withdrawal
        else
            return params["adj_cost_fixed"] + params["adj_cost_prop"] * early_withdrawal
        end 

    else
        return 0.0
    end
end
