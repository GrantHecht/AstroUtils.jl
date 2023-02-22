
struct CubicSpline{LinRangeType <: LinRange}
    # Cubic spline time stamps
    ts::LinRangeType

    # Spline coefficients
    coeffs::Matrix{Float64}

    function CubicSpline(ts, states)
        # Get number of points
        nPoints = length(ts)

        # Allocate matrix of coefficients
        n       = nPoints - 1
        coeffs  = zeros(4*n, 6)

        # Allocate A and b matricies
        A       = spzeros(4*n, 4*n)
        b       = zeros(4*n, 1)

        # Loop through position states and compute splines
        for i in 1:6
            # Fill A and B with 0th order constraints
            A .= 0.0
            b .= 0.0
            for j in 1:n
                x1 = ts[j]
                x2 = ts[j + 1]
                y1 = states[j, i]
                y2 = states[j + 1, i]

                A[1 + 2*(j - 1), 1 + 4*(j - 1)] = x1*x1*x1
                A[1 + 2*(j - 1), 2 + 4*(j - 1)] = x1*x1
                A[1 + 2*(j - 1), 3 + 4*(j - 1)] = x1
                A[1 + 2*(j - 1), 4 + 4*(j - 1)] = 1.0
                A[2 + 2*(j - 1), 1 + 4*(j - 1)] = x2*x2*x2
                A[2 + 2*(j - 1), 2 + 4*(j - 1)] = x2*x2
                A[2 + 2*(j - 1), 3 + 4*(j - 1)] = x2
                A[2 + 2*(j - 1), 4 + 4*(j - 1)] = 1.0
                b[1 + 2*(j - 1)]                = y1
                b[2 + 2*(j - 1)]                = y2
            end

            # Fill A and B with 1st order constraints
            for j in 1:n-1
                x  = ts[j + 1]
                A[2*n + j, 1 + 4*(j - 1)]           = 3.0*x*x
                A[2*n + j, 2 + 4*(j - 1)]           = 2.0*x
                A[2*n + j, 3 + 4*(j - 1)]           = 1.0
                A[2*n + j, 5 + 4*(j - 1)]           = -3.0*x*x
                A[2*n + j, 6 + 4*(j - 1)]           = -2.0*x
                A[2*n + j, 7 + 4*(j - 1)]           = -1.0
            end

            # Fill A and B with 2nd order constriants
            for j in 1:n-1
                x  = ts[j + 1]
                A[2*n + n - 1 + j, 1 + 4*(j - 1)]           = 6.0*x
                A[2*n + n - 1 + j, 2 + 4*(j - 1)]           = 2.0
                A[2*n + n - 1 + j, 5 + 4*(j - 1)]           = -6.0*x
                A[2*n + n - 1 + j, 6 + 4*(j - 1)]           = -2.0
            end

            # If on position state...
            if i < 4
                # Add final constraint to force first and last 1st deriv to equal velocity
                x1  = ts[1]
                x2  = ts[end]
                v1  = states[1, 3 + i]
                v2  = states[end, 3 + i]

                A[end - 1, 1]      = 3.0*x1*x1
                A[end - 1, 2]      = 2.0*x1
                A[end - 1, 3]      = 1.0
                A[end, end - 3]    = 3.0*x2*x2
                A[end, end - 2]    = 2.0*x2
                A[end, end - 1]    = 1.0
                b[end - 1]         = v1
                b[end]             = v2
            else
                # Add final constraint to force first and last 2nd deriv to equal zero
                x1  = ts[1]
                x2  = ts[end]

                A[end - 1, 1]      = 6.0*x1
                A[end - 1, 2]      = 2.0
                A[end, end - 3]    = 6.0*x2
                A[end, end - 2]    = 2.0
            end

            # Set coefficients
            coeffs[:,i] .= A \ b
        end

        return new{typeof(ts)}(ts, coeffs)
    end
end

function getState(spline::CubicSpline, t)
    # Check that t is in [0, 1]
    if t < 0.0 || t > 1.0
        throw(ArgumentError("Time passed to getState() is outside of CubicSpline bounds."))
    end

    # Find relevent polynomial index 
    idxFound = false
    idx      = 0
    while !idxFound
        idx += 1
        if t >= spline.ts[idx] && t <= spline.ts[idx + 1]
            idxFound = true
        end
    end

    # Compute interpolants
    rx  = spline.coeffs[1 + 4*(idx - 1),1]*t^3 + spline.coeffs[2 + 4*(idx - 1),1]*t^2 + spline.coeffs[3 + 4*(idx - 1),1]*t + spline.coeffs[4 + 4*(idx - 1),1]
    ry  = spline.coeffs[1 + 4*(idx - 1),2]*t^3 + spline.coeffs[2 + 4*(idx - 1),2]*t^2 + spline.coeffs[3 + 4*(idx - 1),2]*t + spline.coeffs[4 + 4*(idx - 1),2]
    rz  = spline.coeffs[1 + 4*(idx - 1),3]*t^3 + spline.coeffs[2 + 4*(idx - 1),3]*t^2 + spline.coeffs[3 + 4*(idx - 1),3]*t + spline.coeffs[4 + 4*(idx - 1),3]
    vx  = spline.coeffs[1 + 4*(idx - 1),4]*t^3 + spline.coeffs[2 + 4*(idx - 1),4]*t^2 + spline.coeffs[3 + 4*(idx - 1),4]*t + spline.coeffs[4 + 4*(idx - 1),4]
    vy  = spline.coeffs[1 + 4*(idx - 1),5]*t^3 + spline.coeffs[2 + 4*(idx - 1),5]*t^2 + spline.coeffs[3 + 4*(idx - 1),5]*t + spline.coeffs[4 + 4*(idx - 1),5]
    vz  = spline.coeffs[1 + 4*(idx - 1),6]*t^3 + spline.coeffs[2 + 4*(idx - 1),6]*t^2 + spline.coeffs[3 + 4*(idx - 1),6]*t + spline.coeffs[4 + 4*(idx - 1),6]

    return SVector(rx,ry,rz,vx,vy,vz)
end

function getPosition(spline::CubicSpline, t)
    # Check that t is in [0, 1]
    if t < 0.0 || t > 1.0 # Not working with Symbolics for sparsity detection atm
        throw(ArgumentError("Time passed to getState() is outside of CubicSpline bounds."))
    end

    # Find relevent polynomial index 
    idxFound = false
    idx      = 0
    while !idxFound
        idx += 1
        if t >= spline.ts[idx] && t <= spline.ts[idx + 1]
            idxFound = true
        end
    end

    # Compute interpolants
    rx  = spline.coeffs[1 + 4*(idx - 1),1]*t^3 + spline.coeffs[2 + 4*(idx - 1),1]*t^2 + spline.coeffs[3 + 4*(idx - 1),1]*t + spline.coeffs[4 + 4*(idx - 1),1]
    ry  = spline.coeffs[1 + 4*(idx - 1),2]*t^3 + spline.coeffs[2 + 4*(idx - 1),2]*t^2 + spline.coeffs[3 + 4*(idx - 1),2]*t + spline.coeffs[4 + 4*(idx - 1),2]
    rz  = spline.coeffs[1 + 4*(idx - 1),3]*t^3 + spline.coeffs[2 + 4*(idx - 1),3]*t^2 + spline.coeffs[3 + 4*(idx - 1),3]*t + spline.coeffs[4 + 4*(idx - 1),3]

    return SVector(rx,ry,rz)
end

function getPositionPartial(spline::CubicSpline, t)
    # Check that t is in [0, 1]
    if t < 0.0 || t > 1.0 # Not working with Symbolics for sparsity detection atm
        throw(ArgumentError("Time passed to getState() is outside of CubicSpline bounds."))
    end

    # Find relevent polynomial index 
    idxFound = false
    idx      = 0
    while !idxFound
        idx += 1
        if t >= spline.ts[idx] && t <= spline.ts[idx + 1]
            idxFound = true
        end
    end

    # Compute interpolants
    drxdt  = 3.0*spline.coeffs[1 + 4*(idx - 1),1]*t^2 + 2.0*spline.coeffs[2 + 4*(idx - 1),1]*t + spline.coeffs[3 + 4*(idx - 1),1]
    drydt  = 3.0*spline.coeffs[1 + 4*(idx - 1),2]*t^2 + 2.0*spline.coeffs[2 + 4*(idx - 1),2]*t + spline.coeffs[3 + 4*(idx - 1),2]
    drzdt  = 3.0*spline.coeffs[1 + 4*(idx - 1),3]*t^2 + 2.0*spline.coeffs[2 + 4*(idx - 1),3]*t + spline.coeffs[3 + 4*(idx - 1),3]

    return SVector(drxdt,drydt,drzdt)
end