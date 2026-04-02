# ================================
# Figure 3
# ================================

println("Generating Figure 3...")

# ----------------------
# Environment setup
# ----------------------
using Pkg

project_root = joinpath(@__DIR__, "..")
Pkg.activate(project_root)
Pkg.instantiate()

# ----------------------
# Load packages
# ----------------------
using Plots
using LaTeXStrings
using PolynomialRoots
using NLsolve
using ProgressMeter
using LinearAlgebra
using Statistics
using Measures

# ----------------------
# Paths
# ----------------------
src_path = joinpath(project_root, "src")
results_path = joinpath(project_root, "results")
mkpath(results_path)


# ----------------------
# Load source code
# ----------------------
include(joinpath(src_path, "Extinction.jl"))
include(joinpath(src_path, "DispersionBands.jl"))
include(joinpath(src_path, "FindGammaModes.jl"))
include(joinpath(src_path, "RadiationPatterns.jl"))
include(joinpath(src_path, "parameters.jl"))

# ----------------------
# Radiation pattern 1x1 unit cells (Figure 3(a))
# ----------------------
P1 = RadPat(2,d,1.6,4)
P1 = plot!(subplot=2, xlims=(-1.25*d,1.25*d),left_margin=3mm, right_margin=3mm)
P1 = plot!(title="(a) 1x1 unit cells", titlefontfamily="Computer Modern", subplot=1)
P1 = plot!(title="(b)", titlefontfamily="Computer Modern", subplot=2)

# ----------------------
# Radiation pattern 10x10 unit cells (Figure 3(c))
# ----------------------
beta = 1.6
nl = 20;
np = nl^2;

iB1 = findB1(nl,d,beta)
P2 = RadPat(nl,d,beta,iB1)
P2 = plot!(title="(c) 10x10 unit cells", titlefontfamily="Computer Modern", subplot=1)
P2 = plot!(title="(d)", titlefontfamily="Computer Modern", subplot=2)

# ----------------------
# Radiation pattern 20x20 unit cells (Figure 3(e))
# ----------------------

beta = 1.6
nl = 40;
np = nl^2;

iB1 = findB1(nl,d,beta)
P3 = RadPat(nl,d,beta,iB1)
P3 = plot!(title="(e) 20x20 unit cells", titlefontfamily="Computer Modern", subplot=1)
P3 = plot!(title="(f)", titlefontfamily="Computer Modern", subplot=2)

# ----------------------
# Plot
# ----------------------
P = plot(P1,P2,P3, layout=(3,1), size=(600,1000))

# ----------------------
# Save figure
# ----------------------
savefig(P, joinpath(results_path, "figure3.png"))

