using AstroUtils, StaticArrays, ForwardDiff

# Position and velocity
rv = SVector(6524.834, 6862.875,  6448.296, # [km] 
             4.901327, 5.533756, -1.976341) # [km/s]
μ  = 3.986e5

# Convert to MEE 
mee = convertState(rv, Cartesian, MEE, μ)

# Compute partials with AstroUtils
cartWrtMee = convertStatePartials(mee, MEE, Cartesian, μ)

# Compute partials with ForwardDiff
cartWrtMeeFD = ForwardDiff.jacobian(x -> convertState(x, MEE, Cartesian, μ), mee)

# Compute diff
jacDiff = cartWrtMee .- cartWrtMeeFD
for i in eachindex(jacDiff)
    @test jacDiff[i] ./ eps(cartWrtMee[i]) < 10.0
end