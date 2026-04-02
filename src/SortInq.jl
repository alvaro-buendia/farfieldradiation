
include("eigensolver.jl")
include("GreenMatrix.jl")
include("2DSSHlattice.jl")

# Sorts bulk modes in q-space

function bulkq(nl,d,beta)
np = nl^2;

p0 = pos(nl,d,beta)

nA = findall(x->x==1, [isodd(i)*(iseven(floor((i-1)/nl))) for i in 1:np]); #A sites
nB = findall(x->x==1, [iseven(i)*(iseven(floor((i-1)/nl))) for i in 1:np]); #B sites
nC = findall(x->x==1, [isodd(i)*(isodd(floor((i-1)/nl))) for i in 1:np]); #C sites
nD = findall(x->x==1, [iseven(i)*(isodd(floor((i-1)/nl))) for i in 1:np]); #D sites

evez = eigvecs(Gf(nl,d,beta,wsp0))
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

nbA = nb[findall(x->x in nA, nb)]
nbB = nb[findall(x->x in nB, nb)]
nbC = nb[findall(x->x in nC, nb)]
nbD = nb[findall(x->x in nD, nb)]

# Calculate dispersion bulk bands

 WkA(kx,ky) = [sum(exp.(im*kx*p0[nbA,1]).*exp.(im*ky*p0[nbA,2]).*evez[nbA,n]) for n in ibs] # Fourier transform for sublattice A in bulk sites
 WkB(kx,ky) = [sum(exp.(im*kx*p0[nbB,1]).*exp.(im*ky*p0[nbB,2]).*evez[nbB,n]) for n in ibs] # Fourier transform for sublattice B in bulk sites
 WkC(kx,ky) = [sum(exp.(im*kx*p0[nbC,1]).*exp.(im*ky*p0[nbC,2]).*evez[nbC,n]) for n in ibs] # Fourier transform for sublattice C in bulk sites
 WkD(kx,ky) = [sum(exp.(im*kx*p0[nbD,1]).*exp.(im*ky*p0[nbD,2]).*evez[nbD,n]) for n in ibs] # Fourier transform for sublattice D in bulk sites

# Fourier transforms in the k-path for each sublattice

 WkpA = [WkA(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpB = [WkB(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpC = [WkC(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpD = [WkD(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]

# Rewrite in matrix form

 WkpA = [WkpA[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]
 WkpB = [WkpB[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]
 WkpC = [WkpC[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]
 WkpD = [WkpD[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]

# Classify modes depending on the relative phase between sublattices

Wkp1 = (WkpA+WkpB+ WkpC+WkpD)
Wkp1sq= real(Wkp1.*conj.(Wkp1))

Wkp2 = (WkpA-WkpB+ WkpC-WkpD)
Wkp2sq = real(Wkp2.*conj.(Wkp2))

Wkp3 = (WkpA+WkpB-WkpC-WkpD)
Wkp3sq= real(Wkp3.*conj.(Wkp3))

Wkp4 = (WkpA-WkpB-WkpC+WkpD)
Wkp4sq = real(Wkp4.*conj.(Wkp4))

symw= [findmax([Wkp1sq[i,j]; Wkp2sq[i,j]; Wkp3sq[i,j]; Wkp4sq[i,j]])[2] for i in 1:length(ibs), j in 1:150]; # classify modes depending in the symmetry (1 = s, 2=p_x, 3=p_y, 4=d_xy).

# Just sum over the subset of modes of each symmetry/relative phase

 Wkp1sq= [real(Wkp1sq[i,j]).*(symw[i,j]==1) for i in 1:length(ibs), j in 1:150];
 Wkp2sq = [real(Wkp2sq[i,j]).*(symw[i,j]==2) for i in 1:length(ibs), j in 1:150];
 Wkp3sq = [real(Wkp3sq[i,j]).*(symw[i,j]==3) for i in 1:length(ibs), j in 1:150];
 Wkp4sq = [real(Wkp4sq[i,j]).*(symw[i,j]==4) for i in 1:length(ibs), j in 1:150];

[Wkp1sq, Wkp2sq, Wkp3sq,  Wkp4sq]
end

# Sorts edge states in q-space

function edgeq(nl,d,beta)
np = nl^2;

p0 = pos(nl,d,beta)

# Edge sites in horizontal and verttical edges for different sublattices

nA = findall(x->x==1, [isodd(i)*(iseven(floor((i-1)/nl))) for i in 1:np]); #A sites
nB = findall(x->x==1, [iseven(i)*(iseven(floor((i-1)/nl))) for i in 1:np]); #B sites
nC = findall(x->x==1, [isodd(i)*(isodd(floor((i-1)/nl))) for i in 1:np]); #C sites
nD = findall(x->x==1, [iseven(i)*(isodd(floor((i-1)/nl))) for i in 1:np]); #D sites

evez = eigvecs(Gf(nl,d,beta,wsp0))
nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

nexA = nex[findall(x->x in nA, nex)]
nexB = nex[findall(x->x in nB, nex)]
nexC = nex[findall(x->x in nC, nex)]
nexD = nex[findall(x->x in nD, nex)]

neyA = nex[findall(x->x in nA, ney)]
neyB = nex[findall(x->x in nB, ney)]
neyC = nex[findall(x->x in nC, ney)]
neyD = nex[findall(x->x in nD, ney)]

# Calculate dispersion bulk bands

WkeA(kx,ky) = [sum(exp.(im*kx*p0[nexA,1]).*exp.(im*ky*p0[nexA,2]).*evez[nexA,n]) for n in ies] # Fourier transform for sublattice A in bulk sites
WkeB(kx,ky) = [sum(exp.(im*kx*p0[nexB,1]).*exp.(im*ky*p0[nexB,2]).*evez[nexB,n]) for n in ies] # Fourier transform for sublattice B in bulk sites
WkeC(kx,ky) = [sum(exp.(im*kx*p0[nexC,1]).*exp.(im*ky*p0[nexC,2]).*evez[nexC,n]) for n in ies] # Fourier transform for sublattice C in bulk sites
WkeD(kx,ky) = [sum(exp.(im*kx*p0[nexD,1]).*exp.(im*ky*p0[nexD,2]).*evez[nexD,n]) for n in ies] # Fourier transform for sublattice D in bulk sites

WkpeA = [WkeA(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
WkpeB = [WkeB(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
WkpeC = [WkeC(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
WkpeD = [WkeD(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]

WkpeA = [WkpeA[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]
WkpeB = [WkpeB[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]
WkpeC = [WkpeC[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]
WkpeD = [WkpeD[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]

Wkpe2A = real(WkpeA.*conj.(WkpeA))/maximum(real(WkpeA.*conj.(WkpeA)))
Wkpe2B = real(WkpeB.*conj.(WkpeB))/maximum(real(WkpeB.*conj.(WkpeB)))
Wkpe2C = real(WkpeC.*conj.(WkpeC))/maximum(real(WkpeC.*conj.(WkpeC)))
Wkpe2D= real(WkpeD.*conj.(WkpeD))/maximum(real(WkpeD.*conj.(WkpeD)))

# Fourier transforms in the k-path for each sublattice

 WkpeA = [WkeA(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpeB = [WkeB(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpeC = [WkeC(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpeD = [WkeD(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]

# Rewrite in matrix form

 WkpeA = [WkpeA[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]
 WkpeB = [WkpeB[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]
 WkpeC = [WkpeC[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]
 WkpeD = [WkpeD[i][j] for j in 1:length(ies), i in 1:1:length(kp[:,1])]

# Classify modes depending on the relative phase between sublattices

Wkpe1 = (WkpeA+WkpeB+WkpeC+WkpeD)
Wkpe1sq= real(Wkpe1.*conj.(Wkpe1))

Wkpe2 = (WkpeA+WkpeB-WkpeC-WkpeD)
Wkpe3sq= real(Wkpe2.*conj.(Wkpe2))

Wkpe3 = (WkpeA-WkpeB+WkpeC-WkpeD)
Wkpe2sq = real(Wkpe3.*conj.(Wkpe3))

Wkpe4 = (WkpeA-WkpeB-WkpeC+WkpeD)
Wkpe4sq = real(Wkpe4.*conj.(Wkpe4))

symwe= [findmax([Wkpe1sq[i,j]; Wkpe2sq[i,j]; Wkpe3sq[i,j]; Wkpe4sq[i,j]])[2] for i in 1:length(ies), j in 1:150]; # classify modes depending in the symmetry (1 = s, 2=p_x, 3=p_y, 4=d_xy).

# Just sum over the subset of modes of each symmetry/relative phase

 Wkpe1sq= [real(Wkpe1sq[i,j]).*(symwe[i,j]==1) for i in 1:Int(length(ies)/2), j in 1:150];
 Wkpe2sq = [real(Wkpe2sq[i,j]).*(symwe[i,j]==2) for i in 1:Int(length(ies)/2), j in 1:150];
 Wkpe3sq = [real(Wkpe3sq[i,j]).*(symwe[i,j]==3) for i in Int(length(ies)/2+1):length(ies), j in 1:150];
 Wkpe4sq = [real(Wkpe4sq[i,j]).*(symwe[i,j]==4) for i in  Int(length(ies)/2+1):length(ies), j in 1:150];

[Wkpe1sq, Wkpe2sq, Wkpe3sq,  Wkpe4sq]
end





function bulkq_nl(nl,d,beta,evez)
np = nl^2;

p0 = pos(nl,d,beta)

nA = findall(x->x==1, [isodd(i)*(iseven(floor((i-1)/nl))) for i in 1:np]); #A sites
nB = findall(x->x==1, [iseven(i)*(iseven(floor((i-1)/nl))) for i in 1:np]); #B sites
nC = findall(x->x==1, [isodd(i)*(isodd(floor((i-1)/nl))) for i in 1:np]); #C sites
nD = findall(x->x==1, [iseven(i)*(isodd(floor((i-1)/nl))) for i in 1:np]); #D sites

nc = [1,nl,nl^2-nl+1,nl^2] # corner lattice sites
ne = sort([2:nl-1;nl+1:nl:nl^2-2nl+1;nl^2-nl+2:nl^2-1;2*nl:nl:nl^2-nl]) # edge lattice sites
nex = sort([2:nl-1;nl^2-nl+2:nl^2-1]) # horizontal edge lattice sites
ney = setdiff(ne,nex) # vertical edge lattice sites
nb = setdiff([i for i in 1:np], [ne;nc]) # bulk lattice sites
ics = sort(sortperm(-[norm(evez[nc,i]) for i in 1:np])[1:4]) # indices of corner eigenstates
ies = sort(sortperm(-[norm(evez[ne,i]) for i in 1:np])[1:4(nl-2)]) # indices of edge eigenstates
ibs = setdiff([i for i in 1:np], [ies;ics]) # indices of bulk eigenstates

nbA = nb[findall(x->x in nA, nb)]
nbB = nb[findall(x->x in nB, nb)]
nbC = nb[findall(x->x in nC, nb)]
nbD = nb[findall(x->x in nD, nb)]

# Calculate dispersion bulk bands

 WkA(kx,ky) = [sum(exp.(im*kx*p0[nbA,1]).*exp.(im*ky*p0[nbA,2]).*evez[nbA,n]) for n in ibs] # Fourier transform for sublattice A in bulk sites
 WkB(kx,ky) = [sum(exp.(im*kx*p0[nbB,1]).*exp.(im*ky*p0[nbB,2]).*evez[nbB,n]) for n in ibs] # Fourier transform for sublattice B in bulk sites
 WkC(kx,ky) = [sum(exp.(im*kx*p0[nbC,1]).*exp.(im*ky*p0[nbC,2]).*evez[nbC,n]) for n in ibs] # Fourier transform for sublattice C in bulk sites
 WkD(kx,ky) = [sum(exp.(im*kx*p0[nbD,1]).*exp.(im*ky*p0[nbD,2]).*evez[nbD,n]) for n in ibs] # Fourier transform for sublattice D in bulk sites

# Fourier transforms in the k-path for each sublattice

 WkpA = [WkA(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpB = [WkB(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpC = [WkC(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]
 WkpD = [WkD(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])]

# Rewrite in matrix form

 WkpA = [WkpA[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]
 WkpB = [WkpB[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]
 WkpC = [WkpC[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]
 WkpD = [WkpD[i][j] for j in 1:length(ibs), i in 1:1:length(kp[:,1])]

# Classify modes depending on the relative phase between sublattices

Wkp1 = (WkpA+WkpB+ WkpC+WkpD)
Wkp1sq= real(Wkp1.*conj.(Wkp1))

Wkp2 = (WkpA-WkpB+ WkpC-WkpD)
Wkp2sq = real(Wkp2.*conj.(Wkp2))

Wkp3 = (WkpA+WkpB-WkpC-WkpD)
Wkp3sq= real(Wkp3.*conj.(Wkp3))

Wkp4 = (WkpA-WkpB-WkpC+WkpD)
Wkp4sq = real(Wkp4.*conj.(Wkp4))

symw= [findmax([Wkp1sq[i,j]; Wkp2sq[i,j]; Wkp3sq[i,j]; Wkp4sq[i,j]])[2] for i in 1:length(ibs), j in 1:150]; # classify modes depending in the symmetry (1 = s, 2=p_x, 3=p_y, 4=d_xy).

# Just sum over the subset of modes of each symmetry/relative phase

 Wkp1sq= [real(Wkp1sq[i,j]).*(symw[i,j]==1) for i in 1:length(ibs), j in 1:150];
 Wkp2sq = [real(Wkp2sq[i,j]).*(symw[i,j]==2) for i in 1:length(ibs), j in 1:150];
 Wkp3sq = [real(Wkp3sq[i,j]).*(symw[i,j]==3) for i in 1:length(ibs), j in 1:150];
 Wkp4sq = [real(Wkp4sq[i,j]).*(symw[i,j]==4) for i in 1:length(ibs), j in 1:150];

[Wkp1sq, Wkp2sq, Wkp3sq,  Wkp4sq]
end


function monopartite_nl(nl,d,evez)
np = nl^2;
p0 = pos(nl,d,1)
Wk(kx,ky) = [sum(exp.(im*kx*p0[:,1]).*exp.(im*ky*p0[:,2]).*evez[:,n]) for n in 1:np]
Wkp = [Wk(kp[i,1],kp[i,2]) for i in 1:length(kp[:,1])];
Wkp = [Wkp[i][j] for j in 1:np, i in 1:1:length(kp[:,1])]
Wkpsq = real(Wkp.*conj.(Wkp));
end