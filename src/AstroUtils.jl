module AstroUtils

using LinearAlgebra
using StaticArrays

include("rotations.jl")
include("coordConversions.jl")

export kep2Cart
export cart2Kep
end
