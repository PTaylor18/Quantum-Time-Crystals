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

gate(::GateName"ZZ_couple"; kwargs...) = gate("ZZ"; kwargs...)

rand(Uniform(-π*1.5,-π*0.5)) # disorder in ϕ
rand(Uniform(-π,π)) # disorder in h

# test coupling indices for ZZ gates
coupling_sequence_test = [(1,2),(3,4)]

# create the TC test unitary
test_circuit = [gatelayer("Rx", 4; (θ=π*0.97)),
        [("ZZ_couple", coupling_sequence_test[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence_test)],
        [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:4]]

auto_correlator_test = [[("CZ", (3, 5))], gatelayer("Rx", 4; (θ=π*0.97)),
        [("ZZ_couple", coupling_sequence_test[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence_test)],
        [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:4],
        [("CZ", (3, 5))]]

# pauli Z on all sites
σx2(ψ::MPS) = measure_pauli(ψ, 2, "X")
σz11(ψ::MPS) = measure_pauli(ψ, 11, "Z")
σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]

function coupling_seq(N)
    sequence  = []
    for i=1:1:N-1
      push!(sequence, (i,i+1))
    end
    return sequence
end

coupling_seq_func_test = coupling_seq(4)

function Mz_evolve(N, nsteps)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    Mz_11_vec = Vector{Float64}();
    Mz_auto_correlator_vec = Vector{Float64}();
    coupling_sequence = coupling_seq(N);
    circuit = Vector[];

    auto_correlator = true

    if auto_correlator
        # first layer - can control which qubit is measured for auto-correlator measurement
        layer = [[("CZ", (11, N+1))], gatelayer("Rx", N; (θ=π*0.97)),
                [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)],
                [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
                [("CZ", (11, N+1))]]

        σz_auto_correlation(ψ::MPS) = measure_pauli(ψ, N+1, "X") # x-axis projection of ancilla qubit (google paper)

        obs = Observer([
          "χs" => linkdims,      # bond dimension at each bond
          "χmax" => maxlinkdim,  # maximum bond dimension
          "σᶻ_auto-correlator" => σz_auto_correlation # <ψ|Ẑ(0)Ẑ(t)|ψ> on selected qubit measured with ancilla qubit
          ])

    else
        # first layer
        layer = [gatelayer("Rx", N; (θ=π*0.97)),
                [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)],
                [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]]

        # define the circuit observer
        if N >= 10
            obs = Observer([
              "χs" => linkdims,      # bond dimension at each bond
              "χmax" => maxlinkdim,  # maximum bond dimension
              "σᶻ(11)" => σz11,      # pauli Z on site 11
              "σᶻ" => σz,            # pauli Z on all sites
              ])
        else
            obs = Observer([
              "χs" => linkdims,      # bond dimension at each bond
              "χmax" => maxlinkdim,  # maximum bond dimension
              "σᶻ" => σz,            # pauli Z on all sites
              ])
        end
    end

    for i=1:nsteps
        append!(t_vec, i)
        push!(circuit, layer)
    end

    noise = false
    if noise
        ψ = runcircuit(circuit; (observer!)=obs, noise = ("amplitude_damping", (γ = 0.01,)))
    else
        ψ = runcircuit(circuit; (observer!)=obs)
    end

    if auto_correlator
        Mz_auto_correlator_res = results(obs, "σᶻ_auto-correlator") # results for <ψ|Ẑ(0)Ẑ(t)|ψ> on qubit 11
        for i=1:length(Mz_auto_correlator_res)
            append!(Mz_auto_correlator_vec, mean(Mz_auto_correlator_res[i]))
        end

        return t_vec, Mz_auto_correlator_vec

    else

        Mz_11_res = results(obs, "σᶻ(11)") # results for σᶻ at qubit 11
        for i=1:length(Mz_11_res)
            append!(Mz_11_vec, mean(Mz_11_res[i]))
        end

        Mz_res = results(obs, "σᶻ") # results for σᶻ
        for i=1:length(Mz_res)
            append!(Mz_vec, mean(Mz_res[i]))
        end

        return t_vec, Mz_vec, Mz_11_vec
    end
end

t_vec, Mz_vec, Mz_11_vec = Mz_evolve(20, 250);
plot(t_vec[1:250], Mz_vec[1:250], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",  legend=false)
plot(t_vec[1:250], Mz_11_vec[1:250], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz at qubit 11",  legend=false)


# |+X> state for ancilla qubit prep?
state(::StateName"X+") = [
  1 / sqrt(2)
  1 / sqrt(2)
]

t_vec, Mz_auto_correlator_vec = Mz_evolve(20, 100);
plot(t_vec[1:100], Mz_auto_correlator_vec[1:100], linetype=:steppre, xaxis="Time T", yaxis="<ψ|Ẑ(0)Ẑ(t)|ψ>",  legend=false)

# 2d circuit using SWAP gates?
# approximate unitary circuit with MPO?

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

# add the disorder
# add decay
# benchmark
# x => RxRzRx ?
# finitetemeperature.jl with pastaq but use β = it for exp(-βH) - Dan will look at this

# machine learning with time crystals
# AFM(J>0) lattice and a lattice with hᶻ >> J
# whenever J >> hᶻ => phase 1
# whenever J << hᶻ => phase 2
# therefore classification

# need to know point of phase transition of time crystal so use machine learning
# have a temporal signal - need signal processing => Recurrent neural network

# process evolution of Hamiltonian => for small over rotation, get time crystal
# at which ϵ does the transition happen?
# need to obtain temporal signal


# use yaopastaq - look in tests folder
# look at phase diagram from dtc - supervised learning for phase boundary
# machine learning - frequencies => PCA
# train model on DTC sequences

#
