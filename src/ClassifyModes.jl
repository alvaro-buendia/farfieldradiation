# ============================================
# Eigenmode classification: bulk / edge / corner
# ============================================

# Requires:
# - eigensolver.jl (for Gf, eigvecs, etc.)

include("eigensolver.jl")

using LinearAlgebra

# --------------------------------------------
# Helper: lattice site classification
# --------------------------------------------
function lattice_sites(nl)
    np = nl^2

    # Corner sites
    nc = [1, nl, nl^2 - nl + 1, nl^2]

    # Edge sites
    ne = sort(vcat(
        2:nl-1,                                 # bottom edge
        nl+1:nl:nl^2 - 2nl + 1,                 # left edge
        nl^2 - nl + 2:nl^2 - 1,                 # top edge
        2nl:nl:nl^2 - nl                       # right edge
    ))

    # Horizontal / vertical split (optional)
    nex = sort(vcat(2:nl-1, nl^2 - nl + 2:nl^2 - 1))
    ney = setdiff(ne, nex)

    # Bulk sites
    nb = setdiff(1:np, vcat(ne, nc))

    return nb, ne, nc, nex, ney
end

# --------------------------------------------
# Core classification (given eigenvectors)
# --------------------------------------------
function classify_modes(evez, nl)
    np = nl^2

    nb, ne, nc, _, _ = lattice_sites(nl)

    # Localization weights
    corner_weight = [norm(evez[nc, i]) for i in 1:np]
    edge_weight   = [norm(evez[ne, i]) for i in 1:np]

    # Identify modes by strongest localization
    ics = sort(sortperm(-corner_weight)[1:4])              # corner modes
    ies = sort(sortperm(-edge_weight)[1:4 * (nl - 2)])     # edge modes
    ibs = setdiff(1:np, vcat(ies, ics))                    # bulk modes

    return ibs, ies, ics
end

# --------------------------------------------
# Linearized case (computes eigenvectors internally)
# --------------------------------------------
function Class_lin(nl, d, beta, w)
    # NOTE: assumes Gf(...) depends on wsp0 globally
    evez = eigvecs(Gf(nl, d, beta, wsp0))
    return classify_modes(evez, nl)
end

# --------------------------------------------
# Non-linearized case (eigenvectors provided)
# --------------------------------------------
function Class_nl(nl, d, beta, evez)
    return classify_modes(evez, nl)
end