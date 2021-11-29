using QuantumOptics
using PyPlot

# Parameters
Ω = 0.5
t = [0:0.1:10;]

# Hamiltonian
b = SpinBasis(1//2)
H = Ω*(sigmap(b) ⊗ sigmam(b) + sigmam(b) ⊗ sigmap(b))

ψ0 = spindown(b) ⊗ spinup(b)
tout, ψτ = timeevolution.schroedinger(t, ψ0, H)

# Reduced density matrix
ρ_red = [ptrace(ψ ⊗ dagger(ψ), 1) for ψ=ψτ]
S = [entropy_vn(ρ)/log(2) for ρ=ρ_red]

figure(figsize=(6, 3))
plot(tout, S)

cd(dirname(@__FILE__));
dir = pwd();
mkpath(string(dir,"\\Graphs\\"))
savefig(string(dir,"\\Graphs\\","\\Entaglement_of_2_qubits.png"))
