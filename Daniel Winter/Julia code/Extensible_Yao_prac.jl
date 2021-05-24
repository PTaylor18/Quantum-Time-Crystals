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
using KrylovKit: eigsolve;

bond(n, i) = sum([put(n, i=>σ) * put(n, i+1=>σ) for σ in (X, Y, Z)]);

heisenberg(n) = sum([bond(n, i)
        for i in 1:n-1]);

h = heisenberg(10);
#w, v = eigsolve(mat(h)
#        ,1, :SR, ishermitian=true)

#Listing 7: Hamiltonian evolution is faster with cache
using BenchmarkTools
te = time_evolve(h, 0.1);
te_cache = time_evolve(cache(h), 0.1);


@btime $(rand_state(10)) |> $te;
@btime $(rand_state(10)) |> $te_cache;

h |> decompose

#Listing 8: Circuit simulation is faster without cache
r= rand_state(10);
@btime r |> $(qft(10));

@btime r|> $(cache(qft(10)));

#Listing 9: 10000-layer VQE
using YaoExtensions
n = 10; depth = 10000;

circuit = dispatch!(
        variational_circuit(n, depth),
        :random);

gatecount(circuit)

nparameters(circuit)

h=heisenberg(n);

#This is very slow on the computer even though it is faster then without Yao
#for i = 1:100
#        _, grad = expect'(h, zero_state(n)=>circuit)
#        dispatch!(-, circuit, 1e-3 * grad)
#        println("Step $i, energy = $(expect(h, zero_state(10)=>circuit))")
#end

# This is very very slow
#grad = faithful_grad(h, zero_state(n)=>circuit; nshots=100);

#Listing 10: The eigendecomposition of a QBIR
O=chain(5, put(5,2=>X), put(5,3=>Y))

E, U = YaoBlocks.eigenbasis(O)

#Listing 11: gradient of a maximum mean discrepancy
target_p = normalize!(rand(1<<5));
