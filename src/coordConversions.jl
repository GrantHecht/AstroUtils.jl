
function cartToKep(x::AbstractArray, μ::AbstractFloat)

    # Position and velocity norms
    r = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    v = sqrt(x[4]^2 + x[5]^2 + x[6]^2)

    # Precompute repeated computations
    rinv  = 1.0 / r
    μinv  = 1.0 / μ
    vsq   = v^2
    rDotV = x[1]*x[4] + x[2]*x[5] + x[3]*x[6]
    μRinv = μ*rinv
    twoπ  = 2.0*π

    # Special case tolerance
    scTol = 1.0e-14;

    # Special case flag
    scFlag = 0

    # Angular momentum
    hv = SVector(x[2]*x[6] - x[3]*x[5],
                 x[3]*x[4] - x[1]*x[6],
                 x[1]*x[5] - x[2]*x[4])
    h  = sqrt(hv[1]^2 + hv[2]^2 + hv[3]^2)

    # Normal vector
    nv = SVector(-hv[2], hv[1], 0.0)
    n = sqrt(nv[1]^2 + nv[2]^2)

    # Eccentricity
    term1 = (vsq - μRinv)*μinv
    term2 = rDotV*μinv
    ev = SVector(term1*x[1] - term2*x[4],
                 term1*x[2] - term2*x[5],
                 term1*x[3] - term2*x[6])
    e = sqrt(ev[1]^2 + ev[2]^2 + ev[3]^2)

    # Specific energy
    enrgy = 0.5*vsq - μRinv 

    # Semiparameter
    if e != 1.0
        a = -μ / (2.0*enrgy)
        p = a*(1 - e^2)
    else
        a = Inf
        p = h^2*μinv
    end

    # Inclination
    i = acos(hv[3]/h)

    # Check if elliptical
    if e > scTol

        # RAAN
        Ω = (nv[2] > 0.0) ? acos(nv[1]/n) : 
            twoπ - acos(nv[1]/n)

        # True anomoly
        eDotR = ev[1]*x[1] + ev[2]*x[2] + ev[3]*x[3]
        ν = (rDotV > 0.0) ? acos(eDotR / (r*e)) :
                twoπ - acos(eDotR / (r*e))

        # Argument of periapsis/True longitude of periapsis
        if i > scTol && abs(i - π) > scTol 
            nDotE = nv[1]*ev[1] + nv[2]*ev[2] + nv[3]*ev[3]
            ω = (ev[3] > 0.0) ? acos(nDotE / (n*e)) : 
                    twoπ - acos(nDotE / (n*e))

            return (SVector(a, e, i, Ω, ω, ν), scFlag)

        else # Elliptical equatorial
            scFlat = 1
            ωt = (ev[2] > 0.0) ? acos(e[1] / e) :
                    twoπ - acos(e[1] / e)

            return (SVector(a, e, i, Ω, ωt, ν), scFlag)
        end

    else # Circular
        if i > scTol && abs(i - π) > scTol # inclined

            # RAAN
            Ω = (nv[2] > 0.0) ? acos(nv[1]/n) : 
                twoπ - acos(nv[1]/n)

            # Argument of latitude
            nDotR = nv[1]*x[1] + nv[2]*x[2] + nv[3]*x[3]
            u = (r[3] > 0.0) ? acos(nDotR / (n*r)) :
                    twoπ - acos(nDotR / (n*r))

            scFlag = 2
            return (SVector(a, e, i, Ω, NaN, u), scFlag)

        else # equitorial

            # True longitude
            λt = (x[2] > 0.0) ? acos(r[1] * rinv) :
                    twoπ - acos(r[1] * rinv)

            scFlag = 3
            return (SVector(a, e, i, NaN, NaN, λt), scFlag)
        end
    end
end
    


    