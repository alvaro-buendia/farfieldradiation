# ============================================
# Dispersion bands: bulk and edge
# ============================================

include(joinpath(@__DIR__, "Sortinq.jl"))
include(joinpath(@__DIR__, "eigensolver.jl"))
include(joinpath(@__DIR__, "ClassifyModes.jl"))

using LinearAlgebra

kp = [1/(d)*[[kx for kx in 0:pi/50:pi];[pi for kx in 0:pi/50:pi];[pi-kx for kx in 0:pi/50:pi]] 1/(d)*[[0 for ky in 0:pi/50:pi];[ky for ky in 0:pi/50:pi];[pi-ky for ky in 0:pi/50:pi]]]

# Calculate bulk bands 

function bulkbands(nl,d,beta) 
np = nl^2;
wvec = wlin(nl,d,beta);
evez = eigvecs(Gf(nl,d,beta,wsp0))
sortq = bulkq(nl,d,beta)
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

wB1 =[real(sum(sortq[1][:,i]/sum(sortq[1][:,i]).*wvec[ibs])) for i in 1:1:150];
wB2 =[real(sum(sortq[2][:,i]/sum(sortq[2][:,i]).*wvec[ibs])) for i in 1:1:150];
wB3 =[real(sum(sortq[3][:,i]/sum(sortq[3][:,i]).*wvec[ibs])) for i in 1:1:150];
wB4 =[real(sum(sortq[4][:,i]/sum(sortq[4][:,i]).*wvec[ibs])) for i in 1:1:150];
[wB1  wB2  wB3  wB4]
end

# Calculate edge bands  

function edgebands(nl,d,beta) 
np = nl^2;
wvec = wlin(nl,d,beta);
evez = eigvecs(Gf(nl,d,beta,wsp0))
sortq = edgeq(nl,d,beta)
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

wE1 =[real(sum(sortq[1][:,i]/sum(sortq[1][:,i]).*wvec[ies[1:Int(length(ies)/2)]])) for i in 1:1:150];
wE2 =[real(sum(sortq[2][:,i]/sum(sortq[2][:,i]).*wvec[ies[1:Int(length(ies)/2)]])) for i in 1:1:150];
wE3 =[real(sum(sortq[3][:,i]/sum(sortq[3][:,i]).*wvec[ies[Int(length(ies)/2+1):length(ies)]])) for i in 1:1:150];
wE4 =[real(sum(sortq[4][:,i]/sum(sortq[4][:,i]).*wvec[ies[Int(length(ies)/2+1):length(ies)]])) for i in 1:1:150];
[wE1  wE2  wE3  wE4]
end


#Calculate bulk bands without linearizing (eigenvectors are provided)

function bulkbands_nl(nl,d,beta,wvec,evez) 
np = nl^2;
sortq = bulkq(nl,d,beta)
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

wB1 =[real(sum(sortq[1][:,i]/sum(sortq[1][:,i]).*wvec[ibs])) for i in 1:1:150];
wB2 =[real(sum(sortq[2][:,i]/sum(sortq[2][:,i]).*wvec[ibs])) for i in 1:1:150];
wB3 =[real(sum(sortq[3][:,i]/sum(sortq[3][:,i]).*wvec[ibs])) for i in 1:1:150];
wB4 =[real(sum(sortq[4][:,i]/sum(sortq[4][:,i]).*wvec[ibs])) for i in 1:1:150];
[wB1  wB2  wB3  wB4]
end

# Q-factor bands

function Qbulkbands(nl,d,beta,wvec,evez) 
np = nl^2;
sortq = bulkq(nl,d,beta)
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

Qvec = -real(wvec)./(2*imag(wvec))

wB1 =[real(sum(sortq[1][:,i]/sum(sortq[1][:,i]).*Qvec[ibs])) for i in 1:1:150];
wB2 =[real(sum(sortq[2][:,i]/sum(sortq[2][:,i]).*Qvec[ibs])) for i in 1:1:150];
wB3 =[real(sum(sortq[3][:,i]/sum(sortq[3][:,i]).*Qvec[ibs])) for i in 1:1:150];
wB4 =[real(sum(sortq[4][:,i]/sum(sortq[4][:,i]).*Qvec[ibs])) for i in 1:1:150];
[wB1  wB2  wB3  wB4]
end

#Calculate bulk bands without linearizing (eigenvectors are provided)

function edgebands_nl(nl,d,beta,wvec,evez) 
np = nl^2;
sortq = edgeq(nl,d,beta)
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

wE1 =[real(sum(sortq[1][:,i]/sum(sortq[1][:,i]).*wvec[ies[1:Int(length(ies)/2)]])) for i in 1:1:150];
wE2 =[real(sum(sortq[2][:,i]/sum(sortq[2][:,i]).*wvec[ies[1:Int(length(ies)/2)]])) for i in 1:1:150];
wE3 =[real(sum(sortq[3][:,i]/sum(sortq[3][:,i]).*wvec[ies[Int(length(ies)/2+1):length(ies)]])) for i in 1:1:150];
wE4 =[real(sum(sortq[4][:,i]/sum(sortq[4][:,i]).*wvec[ies[Int(length(ies)/2+1):length(ies)]])) for i in 1:1:150];
[wE1  wE2  wE3  wE4]
end


function Qedgebands(nl,d,beta,wvec,evez) 
np = nl^2;
sortq = edgeq(nl,d,beta)
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates


Qvec = -real(wvec)./(2*imag(wvec))

wE1 =[real(sum(sortq[1][:,i]/sum(sortq[1][:,i]).*Qvec[ies[1:Int(length(ies)/2)]])) for i in 1:1:150];
wE2 =[real(sum(sortq[2][:,i]/sum(sortq[2][:,i]).*Qvec[ies[1:Int(length(ies)/2)]])) for i in 1:1:150];
wE3 =[real(sum(sortq[3][:,i]/sum(sortq[3][:,i]).*Qvec[ies[Int(length(ies)/2+1):length(ies)]])) for i in 1:1:150];
wE4 =[real(sum(sortq[4][:,i]/sum(sortq[4][:,i]).*Qvec[ies[Int(length(ies)/2+1):length(ies)]])) for i in 1:1:150];
[wE1  wE2  wE3  wE4]
end

# Q factor bands monopartite lattice

function monoQbands(nl,d,wvec,evez) 
np = nl^2;
evez = eigvecs(Gf(nl,d,beta,wsp0))
sortq = monopartite_nl(nl,d,evez)
Qvec = -real(wvec)./(2*imag(wvec))

wm =[real(sum(sortq[:,i]/sum(sortq[:,i]).*Qvec)) for i in 1:1:150];
end
