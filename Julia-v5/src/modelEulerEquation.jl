function getmargutility(params::Dict{String,Float64}, cons)

    if params["gamma"] == 1
        margutils = ones(size(cons))./cons
    else
        margutils = cons.^(-params["gamma"])
    end 

end

function getinversemargutility(params::Dict{String,Float64}, margut)

    if params["gamma"] == 1
        invmargut = ones(size(margut))./margut
    else
        invmargut = margut.^(-1/params["gamma"])
    end 

end

function eulerforzero(params::Dict{String,Float64}, itp, A0::Float64, A1::Float64, Y::Float64)
    #-------------------------------------------------------------------------------#
    #This function returns the following quantity:
    #u'(c_t) - b(1+r)u'(c_t+1)
    #This quantity =0 when the Euler equation u'(c_t) = b(1+r)u'(c_t+1)
    #is satified with equality

    #-------------------------------------------------------------------------------#
    # Get marginal utility of consumption tomorrow, given A1

    # this linearly interpolates the grid of points for EdU tomorrow 
    if linearise == 1
        du1AtA1 = getmargutility(params,itp[A1])
    else
        du1AtA1 = itp[A1]
    end

    ## ------------------------------------------------------------------------
    # Check whether tomorrow's (expected) marginal utility negative If so throw an error

    if (du1AtA1 < 0)
       error("approximated marginal utility is negative")
    end

    ## ------------------------------------------------------------------------
    # Get consumption today and the required output
    todaycons = A0 + Y - A1/(1+params["r"])
    euler = getmargutility(params, todaycons) - (params["beta"] * (1+params["r"]) * du1AtA1)
        
end

function eulerforzero_rearrange(A1::Float64, params::Dict{String,Float64}, EdU1_at_A1, A0::Float64, Y::Float64)
    eulerforzero(params, EdU1_at_A1, A0, A1, Y)
end
