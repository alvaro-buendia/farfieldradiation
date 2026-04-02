# ================================
# Figure 8
# ================================

println("Generating Figure 8...")

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

tau = Inf #no losses

# Find eigenmodes and eigenfrequencies without linearizing 
println("Disclaimer: The size of the array is reduced to 5x5 unit cells to allow for fast computation. To reproduce the figure in the paper, set nl = 30 in Figure8.jl")

println("Topological lattice")

nl = 10 # change to nl = 30 to reproduce paper
beta = 1.6

eigsols = eigensolve_nl(nl,d,beta)

wvec = eigsols[1]
evez = eigsols[2]

# Bulk Q-factors
Qb = Qbulkbands(nl,d,beta,wvec,evez) 

# Edge Q-factors
Qe = Qedgebands(nl,d,beta,wvec,evez) 

# Corner Q-factors 
ic = Class_nl(nl,d,beta,evez)[3]
Qvec = -real(wvec)./(2*imag(wvec))

Qc = Qvec[ic]

println("Monopartite lattice")
# Monopartite lattice

eigsols_m = eigensolve_nl(nl,d,1)

wvec_m = eigsols_m[1]
evez_m = eigsols_m[2]

Qm = abs.(monoQbands(nl,d,wvec_m,evez_m))

ll1 = 1+real(ksp)*d/pi*50 #light lines at resonance frequency
ll2 = 150 - real(ksp)*d/pi*50/sqrt(2) #light lines at resonance frequency

P1 = scatter(log10.(Qm), c = :black, markershape=[:circle  :diamond], label=L"\beta = 1", legendfontsize=14, markerstrokewidth=-1, markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
scatter!(log10.(Qb[:,1]), c = :blue, markerstrokewidth=-1,label=L"\beta = 1.6, B1", markershape=:circle, markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
scatter!(log10.(Qb[:,4]), c = :blue, markerstrokewidth=-1,label=L"\beta = 1.6, B4", markershape=:diamond, markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
vline!([ll1,ll2], c=:black, lw=2, ls=:dash,label=L"")

P2 =scatter(log10.(Qm), c = :black, markershape=[:circle  :diamond], label=L"\beta = 1",  legendfontsize=14, markerstrokewidth=-1, markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
scatter!(log10.(Qe[:,1]), c = :lime, markerstrokewidth=-1,  markershape=:circle , label=L"\beta = 1.6, B1", markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
scatter!(log10.(Qe[:,4]), c = :lime, markerstrokewidth=-1,label=L"\beta = 1.6, E4", markershape=:diamond, markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
vline!([ll1,ll2], c=:black, lw=2, ls=:dash, label=L"")

P3 = scatter(log10.(Qm), c = :black, markershape=[:circle  :diamond],markerstrokewidth=-1,legendfontsize=14, label=L"\beta = 1",  markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
scatter!([1:150], log10.(Qc), c = :red, markerstrokewidth=-1, markershape=[:circle :diamond], label=L"\beta = 1.6,\; C1-4", markerlegend=false,tickfontsize=15,tickfontfamily="Computer Modern",framestyle=:box, grid=false, xticks=([1,51,101,151], [L"\Gamma", L"X",L"M",L"\Gamma"]),xlims=(1,151))
vline!([ll1,ll2], c=:black, lw=2, ls=:dash,label=L"")

P = plot(P1,P2,P3, layout=(3,1),size=(600,1200), left_margin=3mm)

# ----------------------
# Save figure
# ----------------------
savefig(P, joinpath(results_path, "figure8.png"))
