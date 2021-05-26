#Listing 6: Heisenberg Hamiltonian
using Yao
using KrylovKit: eigsolve;

bond(n, i) = sum([put(n, i=>σ) * put(n, i+1=>σ) for σ in (X, Y, Z)])

heisenberg(n) = sum([bond(n, i) for i in 1:n-1])

h = heisenberg(10);
w, v = eigsolve(mat(h), 1, :SR, ishermitian=true)
