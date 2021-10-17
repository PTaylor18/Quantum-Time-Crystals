using PastaQ
using ITensors
using Random
using Plots
using Statistics
using Distributions

Random.seed!(1234)
N = 10
B = 1.0

N = 4     # Number of qubits
depth = 10 # Depth of the quantum circuit

# Generate random quantum circuit built out of
# random 2-qubit unitaries
println("Random circuit of depth $depth on $N qubits:")
circuit = randomcircuit(N, depth)

gates = [("X" , 1),                        # Pauli X on qubit 1
         ("CX", (1, 3)),                   # Controlled-X on qubits [1,3]
         ("Rx", 2, (θ=0.5,)),              # Rotation of θ around X
         ("Rn", 3, (θ=0.5, ϕ=0.2, λ=1.2)), # Arbitrary rotation with angles (θ,ϕ,λ)
         ("√SWAP", (3, 4)),                # Sqrt Swap on qubits [2,3]
         ("T" , 4)]                        # T gate on qubit 4

# Build Ising Hamiltonian
sites = siteinds("Qubit", N)
ampo = AutoMPO()
for j in 1:(N - 1)
  # Ising ZZ interactions
  ampo .+= -1, "Z", j, "Z", j + 1
end
for j in 1:N
  # Transverse field X
  ampo .+= -B, "X", j
end
H = MPO(ampo, sites)

# Run DMRG
dmrg_iter = 5      # DMRG steps
dmrg_cutoff = 1E-10   # Cutoff
Ψ0 = randomMPS(sites) # Initial state
sweeps = Sweeps(dmrg_iter)
maxdim!(sweeps, 10, 20, 30, 40, 50, 100)
cutoff!(sweeps, dmrg_cutoff)
# Run DMRG
println("Running DMRG to get ground state of transverse field Ising model:")
E, Ψ = dmrg(H, Ψ0, sweeps)
@show maxlinkdim(Ψ)
println()

nshots = 10_000
data = getsamples(Ψ, nshots; local_basis=["X", "Y", "Z"])

#---

# Define custom function to measure an observable, in this
# case a Pauli operator on `site`
function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

import PastaQ: gate

# Ising (ZZ) coupling gate
function gate(::GateName"ZZ"; ϕ::Number)
  return [
    exp(-im*ϕ/2) 0 0 0
    0 exp(im*ϕ/2) 0 0
    0 0 exp(im*ϕ/2) 0
    0 0 0 exp(-im*ϕ/2)
  ]
end

gate(::GateName"ZZ_couple") = gate("ZZ", ϕ=-0.4)

rand(Uniform(-π*1.5,-π*0.5)) # disorder in ϕ
rand(Uniform(-π,π)) # disorder in h

# coupling indices for ZZ gates
coupling_sequence = [(1,2),(3,4)]

# create the TC Unitary
circuit = [gatelayer("Rx", 4; (θ=π*0.97)),
          [("ZZ_couple", coupling_sequence[i]) for i=1:length(coupling_sequence)]]

# pauli Z on all sites
σx2(ψ::MPS) = measure_pauli(ψ, 2, "X")
σz11(ψ::MPS) = measure_pauli(ψ, 11, "Z")
σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]

function coupling_seq(N)
    sequence  = []
    for i=1:2:N-1
      push!(sequence, (i,i+1))
    end
    return sequence
end

function Mz_evolve(N, nsteps)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    coupling_sequence = coupling_seq(N);
    circuit = Vector[];
    # first layer
    layer = [gatelayer("Rx", N; (θ=π*0.97)),
            [("ZZ_couple", coupling_sequence[i]) for i=1:length(coupling_sequence)],
            gatelayer("Rz", N; (ϕ=π))]
    for i=1:nsteps
        append!(t_vec, i)
        push!(circuit, layer)
    end

    # define the Circuit observer
    obs = Observer([
      "χs" => linkdims,      # bond dimension at each bond
      "χmax" => maxlinkdim,  # maximum bond dimension
      "σᶻ(11)" => σz11,        # pauli X on site 2
      "σᶻ" => σz,
    ])
    ψ = runcircuit(circuit; (observer!)=obs)
    Mz_res = results(obs, "σᶻ")
    for i=1:length(Mz_res)
        append!(Mz_vec, mean(Mz_res[i]))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(50, 100);
plot(t_vec[1:100], Mz_vec[1:100], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",  legend=false)

# 2d circuit using SWAP gates?

function zz_evolve(n::Int, angle::Float64)
    cirq = chain(n)
    for i=1:n-1
        # implement ZZ(θ) unitaries
        push!(cirq, chain(n, cnot(i,i+1)))
        push!(cirq, chain(n, put((i+1)=>Rz(2*angle))))
        push!(cirq, chain(n, cnot(i,i+1)))
    end
    return cirq
end
