function setpar(;beta = 0.95, r_b = 0.04, r = 0.04, adj_cost = true, adj_cost_fixed = 0.1, 
    det_inc = true)

    # Define the parameters as a dictionary
    # TODO: create a Dict of various objects (params, objs, etc)
    # TODO: https://docs.julialang.org/en/v1/manual/performance-tips/
    params                     = Dict{String, Float64}()
    params["tol"]              = 1e-5               # max allowed error
    params["minCons"]          = 1e-5                # min allowed consumption

    # Returns
    # params["r_b"]              = 1.0/0.98 - 1.0      # Interest rate
    params["r_b"]              = r_b
    params["r"]                = r     # Interest rate

    # Preferences 
    params["beta"]             = beta                # 1/(1+r) # Discount factor
    params["gamma"]            = 1.5                 # Coefficient of relative risk aversion
    params["gamma_mod"]        = 1.0-params["gamma"] # For speed, just do this once
    params["startA"]           = 0.0                 # How much asset do people start life with

    # Income process
    params["mu"]               = 0.0                 # mean of initial log income
    params["sigma"]            = 0.2                 # variance of innovations to log income
    # params["sigma"]            = 0.01                # variance of innovations to log income
    params["rho"]              = 0.75                 # persistency of log income
    params["Yretire"]          = 0.5

    # Retirement account
    if adj_cost
        # params["adj_cost_fixed"]   = 0.1                 # fixed cost to adjust the illiquid asset - about 5% of avg annual income as the fixed cost
        # params["adj_cost_prop"]    = 0.1                 # proportional cost to adjust the illiquid asset
        params["adj_cost_fixed"]   = adj_cost_fixed
        params["adj_cost_prop"]    = 0.4
    else
        # No adj costs:
        params["adj_cost_fixed"]   = 0.0 
        params["adj_cost_prop"]    = 0.0
    end

    # Income Polynomial: same as Kovacs Moran
    params["inc_reg_constant"] = 8.200711172
    if det_inc
        params["inc_reg_age"]      = 0.13785127489
        params["inc_reg_age2"]     = -0.00193768738
        params["inc_reg_age3"]     = 6.589127792278049e-6
        # params["Y0_sigma"]         = sqrt( 0.1683 )
    else
        params["inc_reg_age"]      = 0.0
        params["inc_reg_age2"]     = 0.0
        params["inc_reg_age3"]     = 0.0
    end

    return params

        # Eirik Version:
    # function benchpar(;β = 0.98,χ = 0.001,α=1.0,γ=2.0,κ=0.048,d=0.2,
    #     lc=0.0,lb=0.0,ls=0.0,
    #     nx = 41,ns=2,nh=2,ne=1,np=3,nν=3,nϵ=1,
    #         ms = 0.0,mb=0.0,welfare = true,
    #         contract_type = "convex",)

    # ## Initialize housing and ownership grids
    # if ns == 1
    #     sgrd = [0.0] # Only renting
    # elseif ns == 2
    #     sgrd = [0.0, 1.0] # Rent or own
    # else
    #     sgrd = [0.0; 1.0; collect(range(0.5,length=ns-2,step=0.5/(ns-2)))]
    #     # Note that the grid is NOT sorted. 0 (first) = rent, 1 (second) = own, rest is partial ownershi
    # end

    # if nh > 1
    #     hgrd = range(start=1.0,stop=2.0,length=nh)
    # else
    #     hgrd = [1.0]
    # end
    # shgrd = fill(true,(ns,nh))



    # if welfare == true
    #     rentsupport = 10.0 # On top of rent
    #     sellersupport = 10.0 # On top of rent after selling
    # else
    #     rentsupport = -Inf64 # On top of rent
    #     sellersupport = -Inf64 # On top of rent after mincost_selling_to_rent
    # end

    # Eirik approach: named tupple
    # par = (γ=γ,
    #     β = β,
    #     χ = χ,
    #     η = 0.25,
    #     κ = κ,
    #     d = d,
    #     b2y = 5.0,
    #     δ = 0.025,
    #     rf = 0.022,
    #     rm = 0.01,
    #     ρ = 0.95,
    #     α = α,

    #     ms = ms,
    #     mb = mb,

    #     ls = ls, #Fixed cost of selling shared ownership (s∈(0,1))
    #     lc = lc, #Fixed cost of changing ownershare (s∈(0,1)!=s' & h==h')
    #     lb = lb, #Fixed cost of buying shared ownership (s'∈(0,1))

    #     Ts = 25,
    #     Tr = 67,
    #     Te = 100,

    #     xgrd = range(start=-8.0,stop=1000.0,length=nx), # Wealth
    #     hgrd = hgrd, # Contains all house sizes
    #     sgrd = sgrd, # 0=rent, 1=own.
    #     shgrd = shgrd, # Grid for housing×ownership sets (row is tenure (s), column is house quality (h))
        
    #     νgrd = fill(1.0,nν), # Persistent Income Shock (productivity)
        

    #     ## Simulation
    #     N = 6000,
    #     seed = 221956, # Serial Number

    #     ## Various parameters
    #     minrentalsize = minimum(hgrd[shgrd[1,:]]), # Find the smallest rental unit\
    #     rentsupport = rentsupport, # On top of rent
    #     sellersupport = sellersupport, # On top of rent after selling

    #     )


end

# Immutable struct:
# par = benchpar()
# par.β

# Problem: this doesnt work to change values
# par.β = 0.95

