# ============================================
# Far-field radiation patterns
# ============================================

include("2DSSHlattice.jl")
include("eigensolver.jl")

using LinearAlgebra
using Plots

# --------------------------------------------
# Radiation pattern
# --------------------------------------------
"""
    RadPat(nl, d, beta, i)

Compute and plot the far-field radiation pattern
for eigenmode `i`.

Outputs:
- 3D radiation pattern
- Real-space eigenmode distribution
"""
function RadPat(nl, d, beta, i)

    np = nl^2

    # Lattice positions
    p0 = pos(nl, d, beta)

    # Eigenmode
    evec = eigvecs(Gf(nl, d, beta, wsp0))[:, i]

    # ----------------------------------------
    # Far-field structure factor
    # ----------------------------------------
    function Wkqbb(kx, ky, kz)
        phase = exp.(im .* (kx .* p0[:,1] .+
                           ky .* p0[:,2] .+
                           kz .* p0[:,3]))
        return norm(sum(phase .* evec))
    end

    # Wavevector magnitude
    kb = real(k0(wsp0))

    # Spherical parametrization
    kx0(θ, φ) = sin(θ) * cos(φ)
    ky0(θ, φ) = sin(θ) * sin(φ)
    kz0(θ, φ) = cos(θ)

    θ = range(0, π, length=100)
    φ = range(0, 2π, length=100)

    kx = kx0.(θ', φ)
    ky = ky0.(θ', φ)
    kz = kz0.(θ', φ)

    # ----------------------------------------
    # Radiation intensity
    # ----------------------------------------
    Er = norm.(Wkqbb.(kb .* kx, kb .* ky, kb .* kz)) .* (kx.^2 .+ ky.^2)

    # Convert to Cartesian surface
    Erx = Er .* kx
    Ery = Er .* ky
    Erz = Er .* kz

    # ----------------------------------------
    # Normalize axis scaling
    # ----------------------------------------
    xmin, xmax = extrema(Erx)
    ymin, ymax = extrema(Ery)
    zmin, zmax = extrema(Erz)

    Δ = maximum([xmax - xmin, ymax - ymin, zmax - zmin])

    xmid = (xmin + xmax) / 2
    ymid = (ymin + ymax) / 2
    zmid = (zmin + zmax) / 2

    # ----------------------------------------
    # Plot
    # ----------------------------------------
    plot(layout=(1,2))

    # 3D radiation pattern
    surface!(
        subplot=1,
        Erx, Ery, Erz,
        xlims=(xmid - Δ/2, xmid + Δ/2),
        ylims=(ymid - Δ/2, ymid + Δ/2),
        zlims=(zmid - Δ/2, zmid + Δ/2),
        c=cgrad([ :red, :black, :lime],[0,0.5-0.5*zmin/(zmid- Δ/2),0.5,0.5+0.5*zmax/(zmid+Δ/2),1]),
        ticks=false,
        grid=false,
        axis=false,
        aspect_ratio=:equal,
        cb=false
    )

    # Real-space mode profile
    scatter!(
        subplot=2,
        p0[:,1], p0[:,2],
        zcolor=real(evec),
        c=cgrad([:red, :black, :lime]),
        markerstrokewidth=-1,
        markersize=40 / nl,
        clims=(-maximum(abs.(real(evec))), maximum(abs.(real(evec))))
    )

    plot!(
        subplot=2,
        aspect_ratio=:equal,
        legend=false,
        grid=false,
        axis=false,
        cb=false
    )
end