# ============================================
# Eigenvalue solver (linear & nonlinear)
# ============================================

include(joinpath(@__DIR__, "GreenMatrix.jl"))
include(joinpath(@__DIR__, "parameters.jl"))

using LinearAlgebra
using NLsolve
using ProgressMeter
using PolynomialRoots

# --------------------------------------------
# Linearized spectrum
# --------------------------------------------
"""
    wlin(nl, d, beta)

Compute linearized eigenfrequencies using the Green matrix
evaluated at the single-particle resonance frequency.
"""
function wlin(nl, d, beta)
    np = nl^2

    # Green matrix at resonance
    Gsp = Gf(nl, d, beta, wsp0)
    evals = eigvals(Gsp)

    wvec = ComplexF64.(zeros(np))

    # Equation to solve
    eq(w) = 1e-7 .* (1 ./ alph(w) .- evals)

    println("Computing linearized spectrum")

    @showprogress for i in 1:np
        function f!(F, x)
            w = x[1] + im * x[2]
            eq0 = eq(w)[i]
            F[1] = real(eq0)
            F[2] = imag(eq0)
        end

        sol = nlsolve(f!, [wsp0, -1 / tau])
        w = sol.zero[1] + im * sol.zero[2]

        wvec[i] = w / wsp0
    end

    return wvec
end

# --------------------------------------------
# Linear eigenvectors
# --------------------------------------------
"""
    evecs_lin(nl, d, beta)

Compute eigenvectors at resonance (linear approximation).
"""
function evecs_lin(nl, d, beta)
    return eigvecs(Gf(nl, d, beta, wsp0))
end

# --------------------------------------------
# Full nonlinear eigenproblem
# --------------------------------------------
"""
    eigensolve_nl(nl, d, beta)

Solve the full eigenvalue problem:
returns eigenfrequencies and eigenvectors.
"""
function eigensolve_nl(nl, d, beta)
    np = nl^2

    Gzz(w) = Gf(nl, d, beta, w)

    wvec = ComplexF64.(zeros(np))

    # Nonlinear equation
    eq(w) = 1e-7 .* (1 ./ alph(w) .- eigvals(Gzz(w)))

    println("Computing non-linearized spectrum")

    @showprogress for i in 1:np
        function f!(F, x)
            w = x[1] + im * x[2]
            eq0 = eq(w)[i]
            F[1] = real(eq0)
            F[2] = imag(eq0)
        end

        sol = nlsolve(f!, [wsp0, -1e-5 * wsp0])
        w = sol.zero[1] + im * sol.zero[2]

        wvec[i] = w / wsp0
    end

    # Eigenvectors at solved frequencies
    println("Calculating eigenvectors")

    evecs = ComplexF64.(zeros(np, np))

    @showprogress for i in 1:np
        G = Gf(nl, d, beta, wvec[i])
        evecs[:, i] = eigvecs(G)[:, i]
    end

    return wvec, evecs
end