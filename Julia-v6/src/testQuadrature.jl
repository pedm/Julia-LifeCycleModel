# GK Quadrature

using QuadGK

@time quadgk(x -> x^2, 0, 1)[1]

reltol = 0.01
abstol = 0.01
function f(x)
    println("Try $x")
    return exp( -4.0 * abs(x - 0.5) )
end

@time quadgk(f, 0., 1.)
@time quadgk(f, 0., 1.; reltol = reltol, abstol=abstol)
@time quadgk(f, 0., 1.; reltol = 0.1, abstol=0.01)


f(x) = exp( -4.0 * abs(x - 0.5) )
@time quadgk(f, 0, 1)[1]


using StatsFuns

μ = 0.0
σ = 1.0
normpdf(μ, σ, -5.)

normcdf(μ, σ, 0.0)
const normBnd              = 3                   # truncate the normal distrib: ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma

f(x) = normpdf(μ, σ, x)
@time quadgk(f, - normBnd * σ, normBnd * σ)

function f_noisy(x)
    println("x = $x")
    return normpdf(μ, σ, x)
end
@time quadgk(f_noisy, - normBnd * σ, normBnd * σ; reltol = 0.1, abstol=0.1)

xs           = collect(- normBnd * σ:0.01:normBnd * σ)
ys           = map(normcdf, xs)
ys_num       = map(max -> quadgk(f, - normBnd * σ, max)[1], xs)
ys_num_quick = map(max -> quadgk(f, - normBnd * σ, max; reltol = 0.1, abstol=0.1), xs)

@time ys_num_quick = map(max -> quadgk(f, - normBnd * σ, max; reltol = 0.1, abstol=0.1)[1], xs)

using Plots
plotly()
plot(xs, ys, label="true", size=(600, 600))
plot!(xs, ys_num, label = "Integrated")
plot!(xs, ys_num_quick, label = "Integrated (Quick)")
gui()



## Distribution of transitory income shocks
σ = 0.1
xs           = collect(- normBnd * σ:0.001:normBnd * σ)
y_pdf = map( Y_trans -> normpdf(-0.5 * σ, σ, Y_trans), xs)
plot(xs, y_pdf, hover = y_pdf, label="pdf", size=(600, 600), title="Density of Transitory Shocks")



## Gause Hermite
Pkg.add("FastGaussQuadrature")
using FastGaussQuadrature
@time nodes, weights = gausshermite( 10 )

# integrates f(x) = x^2 from -1 to 1
@time dot( weights, nodes.^2 )


## QUESTION: can I integrate over A instead? Then I can just interpolate with respect to A. And I won't have to solve it 10x more
# Example: solve for V as I was doing previously
# Create an interpolated V over A
