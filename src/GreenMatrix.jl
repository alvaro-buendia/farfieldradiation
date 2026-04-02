# Define lattice positions

include("2DSSHlattice.jl")

function Gf(nl,d,beta,w)

p0 = pos(nl,d,beta)
np = nl^2;

Rx = zeros(np,np)
Ry = zeros(np,np)
D = zeros(np,np)
 for i in 1:1:np
     for j in i+1:1:np
dpij = p0[i,:]-p0[j,:]
D[i,j]= norm(dpij)
end
end

D = Symmetric(D)

Dc = D + (Inf)*I(np)

Dm1x = 1 ./ Dc
Dm2x = 1 ./ Dc.^2
Dm3x = 1 ./ Dc.^3

Gm1x = Dm1x # long-range term
Gm2x = im*Dm2x # medium-range term
Gm3x = -Dm3x # short-range term

# Green function matrix z mode
exp.(im*k0(w)*D)/(4*pi).*(Gm1x+1/k0(w)*Gm2x+1/k0(w)^2*Gm3x)
end
