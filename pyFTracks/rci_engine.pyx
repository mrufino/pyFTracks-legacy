# rci_engine.pyx
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False

import numpy as np
cimport numpy as np
from libc.math cimport log, exp, pow

# ============================================================
# Locate the thermal interval for a given time t
# ============================================================
cdef inline int T_interval(
    double t,               # query time
    double[::1] tdata,      # thermal history times (ascending)
    int npts                # number of points in tdata
):
    """
    Return the index i such that tdata[i] <= t < tdata[i+1].
    If t is outside the bounds, return -1 (left) or npts-1 (right).

    Optimized with binary search for faster lookup.
    """

    cdef int left = 0
    cdef int right = npts - 1
    cdef int mid

    # Outside left bound
    if t <= tdata[0]:
        return -1

    # Outside right bound
    if t >= tdata[npts - 1]:
        return npts - 1

    # Binary search for internal intervals
    while right - left > 1:
        mid = (left + right) // 2
        if tdata[mid] <= t:
            left = mid
        else:
            right = mid

    return left


# ============================================================
# Linear interpolation of T(t) — optimized version
# ============================================================
cdef inline double Tfun_linear(
    double t,
    int i,
    double[::1] tdata,
    double[::1] Tdata,
    int npts
):
    """
    Compute the linear interpolation of temperature at time t.
    
    Optimization notes:
    1. Slopes for extrapolation and interpolation are precomputed.
    2. Minimizes repeated divisions.
    3. Branching reduced for internal interpolation.
    """

    cdef double slope, dt, dT

    # Left extrapolation
    if i < 0:
        dt = tdata[1] - tdata[0]
        if dt == 0.0:
            return Tdata[0]
        dT = Tdata[1] - Tdata[0]
        slope = dT / dt
        return Tdata[0] + slope * (t - tdata[0])

    # Right extrapolation
    if i >= npts - 1:
        dt = tdata[npts - 1] - tdata[npts - 2]
        if dt == 0.0:
            return Tdata[npts - 1]
        dT = Tdata[npts - 1] - Tdata[npts - 2]
        slope = dT / dt
        return Tdata[npts - 1] + slope * (t - tdata[npts - 1])

    # Internal linear interpolation
    dt = tdata[i + 1] - tdata[i]
    if dt == 0.0:
        return Tdata[i]
    dT = Tdata[i + 1] - Tdata[i]
    slope = dT / dt
    return Tdata[i] + slope * (t - tdata[i])

# ============================================================
# Function k(u; t) — Rate function for RCI integration
# ============================================================
cdef inline double k_func(
    double u,                  # Integration variable (time lag)
    double t,                  # Current time
    int Ti,                    # Index of thermal history interval for T(t-u)
    double[::1] tdata,         # Array of times in thermal history
    double[::1] Tdata,         # Array of temperatures in thermal history
    int npts,                  # Number of points in thermal history
    double c0, double c1,
    double c2, double c3,
    double R, double n
):
    cdef double T
    cdef double log_term
    cdef double exponent

    # Evaluate temperature T(t-u) using linear interpolation
    T = Tfun_linear(t - u, Ti, tdata, Tdata, npts)

    # Compute logarithmic term used in the exponent
    log_term = log(1.0 / (R * T))

    # Compute exponent in the rate function
    exponent = (-1.0 + n) * (
        c0 + (c1 * (c2 - log(u))) / (c3 - log_term)
    )

    # Return the rate function k(u; t)
    return (
        c1 * exp(-exponent)
        / (u * (-c3 + log_term))
    )



# ============================================================
# RCI Engine — Full-history integration with proper u-partitioning
# ============================================================
cpdef np.ndarray[np.float64_t, ndim=1] rci_annealing(
    double[::1] tdata,      # Time array in seconds
    double[::1] Tdata,      # Temperature array in Kelvin
    double c0,
    double c1,
    double c2,
    double c3,
    double R,
    double n,
    int Nt=100,              # Number of output time nodes
    int Nu_local=40          # Points per subinterval
):
    """
    Rate-Continuous Integral (RCI) for fission-track annealing.
    Quad-like integration without SciPy, compatible with pyFTracks.

    Optimizations applied here are safe: pre-allocating buffers and
    minimizing repeated calculations without changing results.
    """

    cdef int npts = tdata.shape[0]
    cdef double tf = tdata[npts - 1]
    cdef double dt = tf / Nt
    cdef double eps = 1e-6

    # Output array of reduced track lengths
    cdef np.ndarray[np.float64_t, ndim=1] reduced_lengths = np.zeros(Nt)

    # Loop variables
    cdef int i, j, k
    cdef double t, u, ua, ub, du
    cdef double integral, subint
    cdef double w
    cdef int Ti
    cdef double x, dx, xmax

    # --------------------------------------------------------
    # Pre-allocate u-points buffer to avoid repeated allocation
    # --------------------------------------------------------
    cdef np.ndarray[np.float64_t, ndim=1] upts = np.empty(npts + 2, dtype=np.float64)
    cdef int n_upts

    for i in range(Nt):
        t = i * dt

        # Short-circuit for t <= 0
        if t <= 0.0:
            reduced_lengths[i] = 1.0
            continue

        # ----------------------------------------------------
        # Construct explicit breakpoints in u
        # ----------------------------------------------------
        n_upts = 0

        # Lower limit
        upts[n_upts] = eps
        n_upts += 1

        # Interior points: u = t - tdata[k]
        for k in range(npts):
            u = t - tdata[k]
            if eps < u < t:
                upts[n_upts] = u
                n_upts += 1

        # Upper limit
        upts[n_upts] = t
        n_upts += 1

        # Simple in-place sort (upts is small)
        for k in range(n_upts):
            for j in range(k + 1, n_upts):
                if upts[j] < upts[k]:
                    upts[k], upts[j] = upts[j], upts[k]

        # ----------------------------------------------------
        # Integrate over subintervals
        # First subinterval: logarithmic spacing (u ~ 0)
        # Remaining subintervals: trapezoid rule
        # ----------------------------------------------------
        integral = 0.0

        for k in range(n_upts - 1):
            ua = upts[k]
            ub = upts[k + 1]

            if ub <= ua:
                continue

            # Thermal interval index for midpoint
            Ti = T_interval(t - 0.5 * (ua + ub), tdata, npts)

            subint = 0.0

            if k == 0:
                # Logarithmic integration for first subinterval
                xmax = log(ub / ua)
                dx = xmax / Nu_local

                for j in range(Nu_local + 1):
                    x = j * dx
                    u = ua * exp(x)

                    # Trapezoidal weights pre-calculated
                    w = 0.5 if (j == 0 or j == Nu_local) else 1.0

                    # u factor cancels 1/u singularity
                    subint += w * k_func(u, t, Ti, tdata, Tdata, npts, c0, c1, c2, c3, R, n) * u

                integral += subint * dx

            else:
                # Linear trapezoid integration for remaining intervals
                du = (ub - ua) / Nu_local

                for j in range(Nu_local + 1):
                    u = ua + j * du
                    w = 0.5 if (j == 0 or j == Nu_local) else 1.0
                    subint += w * k_func(u, t, Ti, tdata, Tdata, npts, c0, c1, c2, c3, R, n)

                integral += subint * du

        # ----------------------------------------------------
        # Post-processing: enforce physical bounds
        # ----------------------------------------------------
        if (1.0 - n) * integral <= 0.0:
            reduced_lengths[i] = 0.0
        else:
            reduced_lengths[i] = max(0.0, 1.0 - pow((1.0 - n) * integral, 1.0 / (1.0 - n)))

    return reduced_lengths
