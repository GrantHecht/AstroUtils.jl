
function cart2Kep(x::AbstractArray, μ::AbstractFloat)

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
        term  = eDotR / (r*e)
        if abs(term) > 1.0
            term = sign(term)
        end
        ν = (rDotV > 0.0) ? acos(term) : twoπ - acos(term)

        # Argument of periapsis/True longitude of periapsis
        if i > scTol && abs(i - π) > scTol 
            nDotE = nv[1]*ev[1] + nv[2]*ev[2] + nv[3]*ev[3]
            ω = (ev[3] > 0.0) ? acos(nDotE / (n*e)) : 
                    twoπ - acos(nDotE / (n*e))

            return (SVector(a, e, i, Ω, ω, ν), scFlag)

        else # Elliptical equatorial
            scFlag = 1
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
    
function kep2Cart(kep, mu)
    # Special case checks
    if kep[2] == 0.0 && kep[3] == 0.0 # Circular equatorial
        kepc    = SVector(kep[1], kep[2], kep[3], 0.0, 0.0, kep[6])
    elseif kep[2] == 0.0 # Circular inclined
        kepc    = SVector(kep[1], kep[2], kep[3], kep[4], 0.0, kep[6])
    elseif kep[3] == 0.0 # Elliptical equatorial
        kepc    = SVector(kep[1], kep[2], kep[3], 0.0, kep[5], kep[6])
    else # Genaric
         kepc    = SVector(kep[1], kep[2], kep[3], kep[4], kep[5], kep[6])
    end
    
    # Compute semi-parameter
    p       = kepc[1]*(1.0 - kepc[2]^2);
    
    # Compute position and velocity in PQW frame
    term1   = 1.0 / (1.0 + kepc[2]*cos(kepc[6]));
    term2   = sqrt(mu / p);
    rpqw    = @SVector [p*cos(kepc[6])*term1, p*sin(kepc[6])*term1, 0];
    vpqw    = @SVector [-term2*sin(kepc[6]), term2*(kepc[2] + cos(kepc[6])), 0];
    
    # Rotate to inertial reference frame
    r       = rot3Vec(rot1Vec(rot3Vec(rpqw, -kepc[5]), -kepc[3]), -kepc[4]);
    v       = rot3Vec(rot1Vec(rot3Vec(vpqw, -kepc[5]), -kepc[3]), -kepc[4]);
        
    cart = SVector(r[1], r[2], r[3], v[1], v[2], v[3]);
end

    