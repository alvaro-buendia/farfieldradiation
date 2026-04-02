# ============================================
# Extinction cross section map
# ============================================

include(joinpath(@__DIR__, "GreenMatrix.jl"))
include(joinpath(@__DIR__, "2DSSHlattice.jl"))

using LinearAlgebra
using ProgressMeter

# --------------------------------------------
# Extinction map
# --------------------------------------------
"""
    extinction_map(wmin, wmax, wst, nl, d, beta)

Compute extinction cross section as a function of:
- frequency (w)
- k-path (Γ → X → M → Γ)

Returns:
- sigma_map (frequency × k-path)
"""
function extinction_map(wmin, wmax, wst, nl, d, beta)

    np = nl^2

    # Green matrix (linearized at resonance)
    Gz = Gf(nl, d, beta, wsp0)

    # Polarizability matrix inverse
    a_array(w) = inv((1 / alph(w)) * I(np) - Gz)

    # Lattice positions
    p0 = pos(nl, d, beta)

    # ----------------------------------------
    # Incident field polarization
    # ----------------------------------------
    function Ek(w, kx, ky)
        kpar2 = kx^2 + ky^2

        if kpar2 < real(k0(w)^2) && kpar2 != 0
            kz = sqrt(k0(w)^2 - kpar2)
            E = [kz * kx / sqrt(kpar2),
                 kz * ky / sqrt(kpar2),
                 -sqrt(kpar2)]
            return normalize(E)

        elseif kpar2 == 0
            return [1.0, 0.0, 0.0]
        else
            return [0.0, 0.0, 0.0]
        end
    end

    # Incident field at lattice sites
    function Einc(w, kx, ky)
        Ez = Ek(w, kx, ky)[3]
        return Ez .* exp.(im .* (kx .* p0[:, 1] .+ ky .* p0[:, 2]))
    end

    # ----------------------------------------
    # k-path (Γ → X → M → Γ)
    # ----------------------------------------
    kvals = 0:π/50:π

    kx = vcat(kvals, fill(π, length(kvals)), π .- kvals)
    ky = vcat(zeros(length(kvals)), kvals, π .- kvals)

    kp = hcat(kx, ky) ./ d
    nk = size(kp, 1)

    # ----------------------------------------
    # Frequency grid
    # ----------------------------------------
    wpar = collect(wmin:wst:wmax)
    nw = length(wpar)

    sigma_map = zeros(nw, nk)

    # ----------------------------------------
    # Main loop
    # ----------------------------------------
    println("Computing extinction map")

    @showprogress for (iw, w) in enumerate(wpar)

        wphys = w * wsp0
        arw = a_array(wphys)

        for ik in 1:nk
            kx, ky = kp[ik, 1], kp[ik, 2]

            Ein = Einc(wphys, kx, ky)

            sigma_map[iw, ik] =
                abs(4π * k0(wphys) *
                    sum(imag.(arw * Ein .* conj.(Ein))))
        end
    end

    return sigma_map
end