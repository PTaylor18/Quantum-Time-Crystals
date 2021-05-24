using Yao,SymEngine


#Listing 1: quantum Fourier transform
cphase(i,j) = control(i, j=> shift(2π/(2^(i-j+1))));

hcphases(n,i) = chain(n, i==j ?
        put(i=>H) : cphase(j,i) for j in i:n);
qft(n) = chain(hcphases(n,i)
        for i in 1:n)

qft(3)


#Listing 2: apply! and pipe
rand_state(3) |> qft(3);
        #same as apply!(rand_state(3),qft(3))

#Listing 3: inspecting gates
@vars θ

shift(θ) |> mat

control(2,1,2=>shift(θ)) |> mat

#Listing 4: gate decomposition
decompose(x::HGate) = Rz(0.5π)*Rx(0.5π)*Rz(0.5π);
decompose(x::AbstractBlock)=chsubblocks(x, decompose.(subblocks(x)));

qft(3) |> decompose

#Listing 5: inverse QFT
iqft(n)=qft(n)';

iqft(3)

#Listing 6: Heisenberg Hamiltonian
using KrylovKit: eigsolve

bond(n,i) = sum([put(n, i=>σ) * put(n, i+1=>σ) for σ in (X,Y,Z)]);

heisenberg(n) = sum([bond(n,i)
        for i in 1:n-1]);

h=heisenberg(16);
w, v = eigsolve(mat(h)
        ,1, :SR, ishermitian=true)
