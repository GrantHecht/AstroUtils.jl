abstract type OrbitalStateRepresentation end
struct Cartesian <: OrbitalStateRepresentation end
struct Keplerian <: OrbitalStateRepresentation end
struct MEE <: OrbitalStateRepresentation end

function convertState(state::AbstractArray, fromState::Type{OrbitalStateRepresentation}, toState::Type{OrbitalStateRepresentation}, mu)
    throw(ArgumentError("State representation conversion not implemented."))
end

convertState(state::AbstractArray, ::Type{Cartesian}, ::Type{Keplerian}, mu) = cart2Kep(state, mu)
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
            tt    = nDotE / (n*e)
            tt    = (tt >  1.0 ? 1.0 : tt)
            tt    = (tt < -1.0 ? -1.0 : tt)
            ω = (ev[3] > 0.0) ? acos(tt) : (twoπ - acos(tt))

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
    
convertState(state::AbstractArray, ::Type{Keplerian}, ::Type{Cartesian}, mu) = kep2Cart(state, mu)
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

convertState(state::AbstractArray, ::Type{Keplerian}, ::Type{MEE}, mu) = kep2Mee(state)
function kep2Mee(kep)
    # Handle case where kep states are NaN (This may not be correct)
    Ω   = isnan(kep[4]) ? 0.0 : kep[4] 
    ω   = isnan(kep[5]) ? 0.0 : kep[5]

    p   = kep[1]*(1.0 - kep[2]^2)
    f   = kep[2]*cos(ω + Ω)
    g   = kep[2]*sin(ω + Ω)
    h   = tan(0.5*kep[3])*cos(Ω)
    k   = tan(0.5*kep[3])*sin(Ω)
    L   = Ω + ω + kep[6]
    
    return SVector(p,f,g,h,k,L)
end

convertState(state::AbstractArray, ::Type{MEE}, ::Type{Keplerian}, mu) = mee2Kep(state)
function mee2Kep(mee)
    a = mee[1] / (1.0 - mee[2]*mee[2] - mee[3]*mee[3])
    e = sqrt(mee[2]*mee[2] + mee[3]*mee[3])
    i = atan(
        2.0*sqrt(mee[4]*mee[4] + mee[5]*mee[5]),
        1.0 - mee[4]*mee[4] - mee[5]*mee[5],
    )
    Ω = atan(mee[5], mee[4])
    ω = atan(
        mee[3]*mee[4] - mee[2]*mee[5],
        mee[2]*mee[4] + mee[3]*mee[5],
    )
    ν = mee[6] - atan(mee[3],mee[2])
    return SVector(a,e,i,Ω,ω,ν)
end

convertState(state::AbstractArray, ::Type{Cartesian}, ::Type{MEE}, mu) = cart2Mee(state, mu)
function cart2Mee(cart, mu)
    r = SA[cart[1], cart[2], cart[3]]
    v = SA[cart[4], cart[5], cart[6]]

    rdv = dot(r, v)
    rmag = norm(r)
    rhat = r / rmag
    hvec = cross(r, v)
    hmag = norm(hvec)
    hhat = hvec / hmag
    vhat = (rmag*v - rdv*rhat) / hmag

    p   = hmag^2 / mu
    k   =  hhat[1]/(1.0 + hhat[3])
    h   = -hhat[2]/(1.0 + hhat[3])
    kk  = k*k
    hh  = h*h
    s2  = 1.0 + hh + kk
    tkh = 2.0*k*h
    ecc = cross(v,hvec) / mu - rhat

    fhat = SA[
        1.0 - kk + hh,
        tkh,
        -2.0*k,
    ]
    ghat = SA[
        tkh,
        1.0 + kk - hh,
        2.0*h,
    ]
    fhat = fhat / s2
    ghat = ghat / s2
    f    = dot(ecc, fhat)
    g    = dot(ecc, ghat)
    L    = atan(rhat[2]-vhat[1],rhat[1]+vhat[2])

    return SA[p,f,g,h,k,L]
end

convertState(state::AbstractArray, ::Type{MEE}, ::Type{Cartesian}, mu) = mee2Cart(state, mu)
function mee2Cart(mee, mu)
    p = mee[1]
    f = mee[2]
    g = mee[3]
    h = mee[4]
    k = mee[5]
    L = mee[6]

    kk = k*k
    hh = h*h
    tkh = 2.0*k*h
    s2 = 1.0 + hh + kk
    cL = cos(L)
    sL = sin(L)
    w = 1.0 + f*cL + g*sL
    r = p / w
    smp = sqrt(mu/p)
    fhat = SA[
        1.0 - kk + hh,
        tkh,
        -2.0*k,
    ]
    ghat = SA[
        tkh,
        1.0 + kk - hh,
        2.0*h,
    ]
    fhat = fhat / s2
    ghat = ghat / s2
    x = r*cL
    y = r*sL
    xdot = -smp*(g + sL)
    ydot = smp*(f + cL)

    return [x*fhat + y*ghat; xdot*fhat + ydot*ghat]
end
