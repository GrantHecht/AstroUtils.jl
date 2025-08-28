abstract type AbstractConicType end
struct Ellipse <: AbstractConicType end
struct Parabola <: AbstractConicType end
struct Hyperbola <: AbstractConicType end

# An implementation of the algorithm proposed by Bruce A. Conway in "An Improved Algorithm
# Due to Laguerre for the Solution of Kepler's Equation". This is equivalent to the method
# employed within EMTG
function kepler(r0,v0,Δt,mu)
    # Initial preliminaries
    r0mag       = norm(r0)
    v0mag       = norm(v0)

    # Step 1: compute alpha
    alpha   = 2.0/r0mag - v0mag*v0mag/mu
    if alpha > 1e-12
        return kepler(r0,r0mag,v0,v0mag,Δt,mu,alpha,Ellipse())
    elseif alpha < 1e-12
        return kepler(r0,r0mag,v0,v0mag,Δt,mu,alpha,Hyperbola())
    else
        return kepler(r0,r0mag,v0,v0mag,Δt,mu,alpha,Parabola())
    end
end

function kepler(r0,r0mag,v0,v0mag,Δt,mu,alpha,type::AbstractConicType)
    # Preliminaries
    alphatol    = 1e-12
    Xtol        = 1e-12
    max_order   = 10
    max_ipo     = 10 # maximum iterations per order
    r0inv       = 1.0 / r0mag
    sqmu        = sqrt(mu)

    # Step 1: finalize alpha computation
    sigma0  = dot(r0,v0)/sqmu
    sqalpha = sqrt(ifelse(alpha > 0.0,alpha,-alpha))

    # Step 2: compute initial guess
    X_new = 0.0
    if alpha > alphatol
        X_new = alpha * sqmu * Δt
    else
        X_new = 0.1 * sqmu * Δt * r0inv
    end
    r = r0mag

    # Step 3: Perform the Laguerre-Conway-Der iteration
    X = 0.0
    N = 2 # current order
    iters_this_N = 0
    total_iters  = 0
    U0 = 0.0; U1 = 0.0; U2 = 0.0; U3 = 0.0;
    while abs(X - X_new) > Xtol && N < max_order
        # Step 3.1: increment the iteration count
        iters_this_N += 1
        total_iters += 1
        if iters_this_N >= max_ipo
            N += 1
            iters_this_N = 0
        end

        # Step 3.2: update X
        X = X_new

        # Step 3.3: compute U0, U1, U2, U3 for candidate
        U0,U1,U2,U3 = computeUs(alpha,sqalpha,X,type)

        # Step 3.4: compute r and sigma
        r = r0mag * U0 + sigma0 * U1 + U2
        sigma = sigma0 * U0 + (1.0 - alpha * r0mag) * U1

        # Step 3.5: compute F, dF, and ddF
        FX = r0mag * U1 + sigma0 * U2 + U3 - sqmu * Δt
        dFX = r
        ddFX = sigma

        # Step 3.6: Laguerre-Conway or Newton iteration depending on the situation
        sgn = sign(dFX)
        denom = abs((N - 1) * (N - 1) * dFX * dFX - N * (N - 1) * FX * ddFX)
        if denom > 0.0
            dX = N*FX / (dFX + sgn * sqrt(denom))
            X -= dX
        else
            dX = FX / dFX
            X -= dX
        end
    end

    # Step 4: find F, G, Ft, Gt
    F = 1.0 - U2 / r0mag
    G = (r0mag * U1 + sigma0 * U2) / sqmu
    Ft = -sqmu / (r0mag * r) * U1
    Gt = 1.0 - U2 / r

    # Step 5: compute final state
    rf = F*r0 + G*v0
    vf = Ft*r0 + Gt*v0

    return rf, vf
end

function computeUs(alpha,sqalpha,X,::Ellipse)
    y = alpha * X * X
    C = (1.0 - cos(sqrt(y))) / y
    S = (sqrt(y) - sin(sqrt(y))) / sqrt(y*y*y)
    U1 = X * (1.0 - y * S)
    U2 = X * X * C
    U3 = X * X * X * S
    U0 = 1.0 - alpha * U2
    return U0, U1, U2, U3
end

function computeUs(alpha,sqalpha,X,::Hyperbola)
    sqalphaX = sqalpha * X
    if abs(sqalphaX) > 30.0
        error("kepler solver failed, went to far out on a hyperbola.")
    end
    U0 = cosh(sqalphaX)
    U1 = sinh(sqalphaX) / sqalphaX
    U2 = (1.0 - U0) / alpha
    U3 = 1.0 / alpha * (X - U1)
    return U0, U1, U2, U3
end

function computeUs(alpha,sqalpha,X,::Parabola)
    U0 = 1.0
    U1 = X
    U2 = 0.5*U1*X
    U3 = U2*X / 3.0
    return U0, U1, U2, U3
end
