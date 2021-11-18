using PastaQ
using ITensors
using Plots
using Statistics
using Distributions
using FFTW

X_flip(N::Int) = gatelayer("Rx", 4; (θ=π*0.97))

# Arrange ZZ(:rand) couplings using CNOT-Rz-CNOT sequence
function ZZ_layer(N::Int)
    qarray = PastaQ.lineararray(N)
    first_sublattice_CNOTs = gatelayer("CX", qarray[1])
    second_sublattice_CNOTs = gatelayer("CX", qarray[2])
    randϕ = rand(Uniform(-π*1.5,-π*0.5), N) # set disordered angles for Ising coupling
    first_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 1:2:N]
    second_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 2:2:N]
    circuit = [[first_sublattice_CNOTs], [first_sublattice_Rzs], [first_sublattice_CNOTs], [second_sublattice_CNOTs], [second_sublattice_Rzs], [second_sublattice_CNOTs]]
    return circuit
end

Z_field(N::Int) = [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]

function entropy_von_neumann(ψ::MPS, b::Int)
  orthogonalize!(ψ, b)
  _, S = svd(ψ[b], (linkind(ψ, b-1), siteind(ψ, b)))
  SvN = 0.0
  for n in ITensors.dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
  end
  return SvN
end

N = 8
s = siteinds("Qubit", N)
psi = randomMPS(s, 4)
SvN = entropy_von_neumann(psi, 4)

linkdims(psi)
maxlinkdim(psi)

psi = randomstate(s, 4; χ = 10)

randomstate
productstate

function entanglemententropy_evolve(N, nsteps, b)
    t_vec = Vector{Float64}();
    SvN_vec = Vector{Float64}();

    #ψ = productstate(N)
    #ψ = randomstate(N)

    s = siteinds("Qubit", N)
    ψ = randomMPS(s, 2)

    for i=1:nsteps
        println("Step: $i")
        append!(t_vec, i)

        ψ = runcircuit(ψ, X_flip(N))
        normalize!(ψ)
        ψ = runcircuit(ψ, ZZ_layer(N))
        normalize!(ψ)
        ψ = runcircuit(ψ, Z_field(N))
        normalize!(ψ)

        SvN = entropy_von_neumann(ψ, b)
        append!(SvN_vec, SvN)
    end
    return t_vec, SvN_vec
end

N = 8 # no. of qubits
nsteps = 200
b = 4 # select site to partition system

t_vec, SvN_vec = entanglemententropy_evolve(N, nsteps, b);
plot(t_vec[1:nsteps], SvN_vec[1:nsteps], xaxis="Time T", yaxis="von Neumann entanglement entropy, SvN", legend=false)

SvN_vec

s = siteinds("Qubit", N)
ψ = productstate(s)


ψ = runcircuit(ψ, X_flip(N))
normalize!(ψ)
ψ = runcircuit(ψ, ZZ_layer(N))
normalize!(ψ)
ψ = runcircuit(ψ, Z_field(N))
normalize!(ψ)

SvN = entropy_von_neumann(ψ, 4)
