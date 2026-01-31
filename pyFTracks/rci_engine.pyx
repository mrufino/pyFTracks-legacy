# rci_engine.pyx
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False

import numpy as np
cimport numpy as np
from libc.math cimport log, exp, pow

# ============================================================
# Localiza o intervalo de interpolação térmica
# ============================================================
cdef inline int T_interval(
    double t,
    double[::1] tdata,
    int npts
):
    cdef int i

    if t <= tdata[0]:
        return -1

    if t >= tdata[npts - 1]:
        return npts - 1

    for i in range(npts - 1):
        if tdata[i] <= t < tdata[i + 1]:
            return i

    return npts - 1


# ============================================================
# Interpolação linear T(t) — versão final
# ============================================================
cdef inline double Tfun_linear(
    double t,
    int i,
    double[::1] tdata,
    double[::1] Tdata,
    int npts
):
    cdef double slope

    # Extrapolação à esquerda
    if i < 0:
        slope = (Tdata[1] - Tdata[0]) / (tdata[1] - tdata[0])
        return Tdata[0] + slope * (t - tdata[0])

    # Extrapolação à direita
    if i >= npts - 1:
        slope = (Tdata[npts - 1] - Tdata[npts - 2]) / (
            tdata[npts - 1] - tdata[npts - 2]
        )
        return Tdata[npts - 1] + slope * (t - tdata[npts - 1])

    # Interpolação linear interna
    return (
        Tdata[i]
        + (Tdata[i + 1] - Tdata[i])
        * (t - tdata[i])
        / (tdata[i + 1] - tdata[i])
    )

# ============================================================
# Função k(u; t) — versão final
# ============================================================
cdef inline double k_func(
    double u,
    double t,
    int Ti,                     # índice do intervalo de T(t-u)
    double[::1] tdata,
    double[::1] Tdata,
    int npts,
    double c0,
    double c1,
    double c2,
    double c3,
    double R,
    double n
):
    cdef double T
    cdef double log_term
    cdef double exponent

    # Avaliação de T(t-u) sabendo o intervalo
    T = Tfun_linear(t - u, Ti, tdata, Tdata, npts)

    log_term = log(1.0 / (R * T))

    exponent = (-1.0 + n) * (
        c0 + (c1 * (c2 - log(u))) / (c3 - log_term)
    )

    return (
        c1 * exp(-exponent)
        / (u * (-c3 + log_term))
    )


# ============================================================
# Motor RCI — versão FINAL com partição correta em u
# ============================================================
cpdef np.ndarray[np.float64_t, ndim=1] rci_annealing(
    double[::1] tdata,      # tempo (s)
    double[::1] Tdata,      # temperatura (K)
    double c0,
    double c1,
    double c2,
    double c3,
    double R,
    double n,
    int Nt=100,
    int Nu_local=40         # pontos por subintervalo
):
    """
    Integral da Constante de Taxa (RCI)
    Versão correta (Mathematica/quad-like), sem SciPy.
    """

    cdef int npts = tdata.shape[0]
    cdef double tf = tdata[npts - 1]
    cdef double dt = tf / Nt
    cdef double eps = 1e-6

    cdef np.ndarray[np.float64_t, ndim=1] reduced_lengths = np.zeros(Nt)

    cdef int i, j, k
    cdef double t, u, ua, ub, du
    cdef double integral, subint
    cdef double w
    cdef int Ti
    cdef double x, dx, xmax


    # buffer fixo para pontos de quebra (npts + eps + t)
    cdef np.ndarray[np.float64_t, ndim=1] upts = np.empty(npts + 2, dtype=np.float64)
    cdef int n_upts

    for i in range(Nt):
        t = i * dt

        if t <= 0.0:
            reduced_lengths[i] = 1.0
            continue

        # ----------------------------------------------------
        # Construção explícita dos pontos de quebra em u
        # ----------------------------------------------------
        n_upts = 0

        # limite inferior
        upts[n_upts] = eps
        n_upts += 1

        # pontos u = t - tdata[k]
        for k in range(npts):
            u = t - tdata[k]
            if eps < u < t:
                upts[n_upts] = u
                n_upts += 1

        # limite superior
        upts[n_upts] = t
        n_upts += 1

        # ordenação simples (n_upts é pequeno e fixo)
        for k in range(n_upts):
            for j in range(k + 1, n_upts):
                if upts[j] < upts[k]:
                    upts[k], upts[j] = upts[j], upts[k]

        # ----------------------------------------------------
        # Integração por subintervalos suaves
        # (primeiro em log(u), restantes em u)
        # ----------------------------------------------------
        integral = 0.0

        for k in range(n_upts - 1):
            ua = upts[k]
            ub = upts[k + 1]

            if ub <= ua:
                continue

            # índice do intervalo térmico (fixo neste subintervalo)
            Ti = T_interval(t - 0.5 * (ua + ub), tdata, npts)

            subint = 0.0

            # -----------------------------
            # PRIMEIRO subintervalo: u ~ 0
            # integração logarítmica
            # -----------------------------
            if k == 0:

                xmax = log(ub / ua)
                dx = xmax / Nu_local

                for j in range(Nu_local + 1):
                    x = j * dx
                    u = ua * exp(x)
                    w = 0.5 if (j == 0 or j == Nu_local) else 1.0

                    # fator u cancela a singularidade 1/u
                    subint += w * k_func(
                        u, t, Ti,
                        tdata, Tdata, npts,
                        c0, c1, c2, c3, R, n
                    ) * u

                integral += subint * dx

            # -----------------------------
            # DEMAIS subintervalos: trapézio normal
            # -----------------------------
            else:
                du = (ub - ua) / Nu_local

                for j in range(Nu_local + 1):
                    u = ua + j * du
                    w = 0.5 if (j == 0 or j == Nu_local) else 1.0

                    subint += w * k_func(
                        u, t, Ti,
                        tdata, Tdata, npts,
                        c0, c1, c2, c3, R, n
                    )

                integral += subint * du

        # ----------------------------------------------------
        # Pós-processamento físico
        # ----------------------------------------------------
        if (1.0 - n) * integral <= 0.0:
            reduced_lengths[i] = 0.0
        else:
            reduced_lengths[i] = max(
                0.0,
                1.0 - pow((1.0 - n) * integral, 1.0 / (1.0 - n))
            )

    return reduced_lengths
