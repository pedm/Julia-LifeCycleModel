using Statistics

"""

Convenience interface for discrete distributions .

If you want to do S{T} <: FiniteSupportDistribution{T}, you need to implement:

- `supp(s::S)` returning a Vector{T} with the support points

  This will fall back to doing `s.support` if you don't define a separate
  method.

- `pmf(s::S, i)` returning the mass of the i-th support element

"""
abstract type FiniteSupportDistribution{T} <: AbstractArray{T, 1} end


Base.size(fsd::FiniteSupportDistribution)=(length(fsd), )
Base.IndexStyle(::Type{<:FiniteSupportDistribution}) = IndexLinear()
Base.getindex(fsd::FiniteSupportDistribution, i)=supp(fsd)[i]



supp(fsd::FiniteSupportDistribution)=fsd.support
Base.length(fsd::FiniteSupportDistribution)=length(supp(fsd))

" Mean of the finite supp. distribution (not the underlying generator) "
Statistics.mean(fsd::FiniteSupportDistribution)=begin
    m = 0.
    for i in 1:length(fsd)
        m += fsd[i] * pmf(fsd, i)
    end
    return m
end



" Divide support points by `mean(fsd)`"
normalized(fsd::FiniteSupportDistribution)=(fsd2 = deepcopy(fsd); normalized!(fsd2); fsd2)
normalized!(fsd::FiniteSupportDistribution)=(rmul!(supp(fsd), 1/mean(fsd)))




struct BasicFiniteSupportDistribution{T} <: FiniteSupportDistribution{T}
    support::Vector{T}
    probvec::Vector{Float64}
end

pmf(fsd::BasicFiniteSupportDistribution)=fsd.probvec
pmf(fsd::BasicFiniteSupportDistribution, i::Int64)=fsd.probvec[i]



struct EquiprobDistribution{T} <: FiniteSupportDistribution{T}
    support::Vector{T}
    equiprob::Float64 # Avoids (super minor) cost of computing 1/n every time
end
EquiprobDistribution(support)=EquiprobDistribution([x for x in support],
                                                   1/length(support))


pmf(epd::EquiprobDistribution)=fill(epd.equiprob,length(epd))
pmf(epd::EquiprobDistribution, i::Int64)=epd.equiprob
Statistics.mean(epd::EquiprobDistribution)=sum(supp(epd))*pmf(epd,1)



struct Discretized1dDistribution{D, S, T <: FiniteSupportDistribution{S}} <: FiniteSupportDistribution{S}
    fsd::T
    distrib::D
end

supp(dd::Discretized1dDistribution)=supp(dd.fsd)
pmf(dd::Discretized1dDistribution)=pmf(dd.fsd)
pmf(dd::Discretized1dDistribution, i)=pmf(dd.fsd,i)
recover(dd::Discretized1dDistribution)=dd.distrib





function eqprob_discretize(dd, npts)
    # Approximate distribution `dd` with a finite support distribution with
    # equally probable points.
    # --------------------------------------------------------------------------
    #
    # Idea of aproximation -----------------------------------------------------
    #     ^
    #     |
    #     |
    # 1   -                                                  o
    #     |
    #     x                                       o
    #     |
    #     -                               o
    #     |
    #     x                           o
    #     |
    #     -                      o
    #     |
    #     x              o
    #     |
    #     -        o
    #     |
    #     x  o
    # 0 _ |__________________________________________________________\
    #                                                                /
    #  Assign F(-) - F(prev(-)) to `x`s
    # --------------------------------------------------------------------------

    ran1 = range(0, (npts-1)/npts, length=npts)
    ran2 = range(1/npts, 1, length=npts)
    ran = (ran1 + ran2)/2

    zsupport = quantile.(dd, ran)

    Discretized1dDistribution(EquiprobDistribution(zsupport), dd)
end
