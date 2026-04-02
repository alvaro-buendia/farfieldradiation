# ============================================
# Gamma-point mode identification
# ============================================

include(joinpath(@__DIR__, "ClassifyModes.jl"))

using LinearAlgebra

# --------------------------------------------
# Sublattice decomposition (A, B, C, D)
# --------------------------------------------
function sublattices(nl)
    np = nl^2

    nA = [i for i in 1:np if isodd(i)  && iseven(fld(i - 1, nl))]
    nB = [i for i in 1:np if iseven(i) && iseven(fld(i - 1, nl))]
    nC = [i for i in 1:np if isodd(i)  && isodd(fld(i - 1, nl))]
    nD = [i for i in 1:np if iseven(i) && isodd(fld(i - 1, nl))]

    return nA, nB, nC, nD
end

# --------------------------------------------
# Symmetry patterns
# --------------------------------------------
function symmetry_vectors(nl)
    np = nl^2
    nA, nB, nC, nD = sublattices(nl)

    sgn1 = [ (i in nA || i in nD ?  1 : 0) - (i in nB || i in nC ? 1 : 0) for i in 1:np ]
    sgn2 = [ (i in nA || i in nC ?  1 : 0) - (i in nB || i in nD ? 1 : 0) for i in 1:np ]
    sgn3 = [ (i in nA || i in nB ?  1 : 0) - (i in nC || i in nD ? 1 : 0) for i in 1:np ]
    sgn4 = [ (i in nA || i in nD ?  1 : 0) + (i in nB || i in nC ? 1 : 0) for i in 1:np ]

    return (sgn1, sgn2, sgn3, sgn4)
end

# --------------------------------------------
# Core selector
# --------------------------------------------
"""
    find_mode(nl, d, beta, indices, symmetry)

Find the mode index within `indices` that maximizes
overlap with a given symmetry vector.
"""
function find_mode(nl, d, beta, indices, symmetry)
    np = nl^2;
    if nl < 2
        return 1
    end

    evecs = eigvecs(Gf(nl, d, beta, wsp0))
    overlaps = abs.(sum(evecs[:, indices] .* symmetry, dims=1))[1,:]
    idx = argmax(overlaps)

    return indices[idx]
end

# --------------------------------------------
# Bulk modes (B1–B4)
# --------------------------------------------
function findB(nl, d, beta, n)
    np = nl^2;
    ibs = Class_lin(nl, d, beta, wsp0)[1]
    sgn = symmetry_vectors(nl)
    evecs = eigvecs(Gf(nl, d, beta, wsp0))
    overlaps = abs.(sum(evecs[:, ibs] .* sgn[n], dims=1))[1,:]
    idx = findmax(overlaps)[2]
    return(ibs[idx])
end

findB1(nl, d, beta) = findB(nl, d, beta, 1)
findB2(nl, d, beta) = findB(nl, d, beta, 2)
findB3(nl, d, beta) = findB(nl, d, beta, 3)
findB4(nl, d, beta) = findB(nl, d, beta, 4)

# --------------------------------------------
# Edge modes (E1–E4)
# --------------------------------------------
function findE(nl, d, beta, n)
    sgn = symmetry_vectors(nl)
    ies = Class_lin(nl, d, beta, wsp0)[2]
    return find_mode(nl, d, beta, ies, sgn[n])
end

findE1(nl, d, beta) = findE(nl, d, beta, 1)
findE2(nl, d, beta) = findE(nl, d, beta, 2)
findE3(nl, d, beta) = findE(nl, d, beta, 3)
findE4(nl, d, beta) = findE(nl, d, beta, 4)

# --------------------------------------------
# Corner modes (C1–C4)
# --------------------------------------------
function findC(nl, d, beta, n)
    sgn = symmetry_vectors(nl)
    ics = Class_lin(nl, d, beta, wsp0)[3]
    return find_mode(nl, d, beta, ics, sgn[n])
end

findC1(nl, d, beta) = findC(nl, d, beta, 1)
findC2(nl, d, beta) = findC(nl, d, beta, 2)
findC3(nl, d, beta) = findC(nl, d, beta, 3)
findC4(nl, d, beta) = findC(nl, d, beta, 4)