#Define positions of 2D SSH lattice, depending on size, period and beta parameters

function pos(nl,d,beta)
np = nl^2;
p0 = zeros(np,3)
    for i in 1:1:np
p0[i,:]= [beta*d/2*floor((mod(i,nl)+nl*iszero(mod(i,nl)))/2)+(d-beta*d/2)*floor(mod(i-1,nl)/2), (beta*d/2)*floor(ceil(i/nl)/2)+(d-beta*d/2)*floor((ceil(i/nl)-1)/2),0]
end
p0[:,1] = p0[:,1] .- mean(p0[:,1]);
p0[:,2] = p0[:,2] .- mean(p0[:,2]);
p0
end
