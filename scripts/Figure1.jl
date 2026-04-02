# ================================
# Figure 1(c)
# ================================

println("Generating Figure 1(c)...")

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
include(joinpath(src_path, "parameters.jl"))
include(joinpath(src_path, "DispersionBands.jl"))

# ----------------------
# Parameters
# ----------------------
wmin = 0.96
wmax = 1.06
wst  = 0.0005

# (must be defined in src or globally)
# d, c, eb, wsp, wsp0, ksp

# ----------------------
# k-path (Γ → X → M → Γ)
# ----------------------
kvals = 0:π/50:π

kx = vcat(kvals, fill(π, length(kvals)), π .- kvals)
ky = vcat(zeros(length(kvals)), kvals, π .- kvals)

kp = hcat(kx, ky) ./ d

# Light line
ll = c * sqrt.(kp[:,1].^2 .+ kp[:,2].^2) ./ sqrt(eb) ./ real(wsp)

# ----------------------
# Extinction map
# ----------------------
sigma_map = extinction_map(wmin, wmax, wst, 30, d, 1.6)

p = heatmap(
    1:length(kx),
    wmin:wst:wmax,
    log10.(real(sigma_map .+ 1e-6)),
    clims = (0, 4.25),
    framestyle = :box,
    grid = false,
    xticks = ([1, 51, 101, 151], [L"\Gamma", L"X", L"M", L"\Gamma"]),
    yticks = wmin:0.02:wmax,
    tickfontsize = 15,
    tickfontfamily = "Computer Modern",
    xlims = (1, 151),
    ylims = (wmin, wmax),
    title = L"\log_{10}(\sigma_\textrm{ext})"
)

# Light cone
plot!(p, ll, ls = :dash, c = :white, lw = 2)

# ----------------------
# Bulk bands
# ----------------------
wbulk = bulkbands(30, d, 1.6)
scatter!(p, 1:5:150, wbulk[1:5:end, :], c = :blue, markerstrokewidth = -1)

# ----------------------
# Edge bands
# ----------------------
wedge = edgebands(30, d, 1.6)
scatter!(p, 1:5:150, wedge[1:5:end, :], c = :lime, markerstrokewidth = -1)

# ----------------------
# Corner modes
# ----------------------
wcorner = wlin(12, d, 1.6)[Class_lin(12, d, 1.6, wsp0)[3]]
scatter!(p, 1:5:150, real(wcorner), c = :red, markerstrokewidth = -1)

# ----------------------
# Annotations
# ----------------------
annotate!(p, 76, 1.047, text(L"B_{4}", :white))
annotate!(p, 76, 0.9775, text(L"B_{1}", :white))
annotate!(p, 76, 0.989, text(L"E_{1,2}", :white))
annotate!(p, 76, 1.022, text(L"E_{3,4}", :white))
annotate!(p, 76, 1.005, text(L"C_{1-4}", :white))
annotate!(p, 101, 0.9975, text(L"B_{2,3}", :white))

# ----------------------
# Inset: Brillouin zone
# ----------------------
x = [-1, 1, 1, -1, -1]
y = [-1, -1, 1, 1, -1]

xk = [0, 1, 1, 0]
yk = [0, 0, 1, 0]

rll = real(ksp) * d / π
θ = range(0, 2π, length = 101)

xcirc = rll .* cos.(θ)
ycirc = rll .* sin.(θ)

plot!(p, inset = (1, bbox(0.15, -0.18, 0.25, 0.25, :center)))

plot!(
    x, y,
    subplot = 2,
    aspect_ratio = :equal,
    xlims = (-1.5, 1.5),
    ylims = (-1.5, 1.5),
    legend = false,
    lw = 1,
    grid = false,
    axis = false,
    ticks = false,
    background_color = :transparent,
    c = :gray
)

plot!(subplot = 2, xcirc, ycirc, fillrange = rll, c = :white)
plot!(subplot = 2, xk, yk, c = :red, lw = 2)
scatter!(subplot = 2, xk, yk, c = :red, markersize = 4, markerstrokewidth = -1)

annotate!(subplot = 2, 1.45, 0.5, text(L"M", :red, 14))
annotate!(subplot = 2, -0.3, 0, text(L"\Gamma", :red, 14))
annotate!(subplot = 2, 1.45, -0.5, text(L"X", :red, 14))
plot!(legend=false)

# ----------------------
# Save figure
# ----------------------
savefig(p, joinpath(results_path, "figure1c.png"))

println("Figure 1(c) saved.")