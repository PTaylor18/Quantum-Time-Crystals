using PastaQ
using ITensors
using Random
using Plots
using Statistics
using Distributions

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

auto_correlator_test = [[("H", 5)],[("CZ", (3, 5))], gatelayer("Rx", 4; (θ=π*0.97)),
        [("ZZ_couple", coupling_sequence_test[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence_test)],
        [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:4],
        [("CZ", (3, 5))],
        [("H", 5)]]

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

# create cat state: |000> + |111>
function cat_state(N)

    circuit = Vector[]
    hadamard = [("H", 1)]
    push!(circuit, hadamard)

    for i=1:N-1
        cnot_layer = [("CNOT",(i,i+1))]
        push!(circuit, cnot_layer)
    end
    return circuit
end

cat_state(4)

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
        layer = [[("H", N+1)],[("CZ", (5, N+1))], gatelayer("Rx", N; (θ=π*0.97)),
                [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)],
                [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
                [("CZ", (5, N+1))],
                [("H", N+1)]]

        σz_auto_correlation(ψ::MPS) = measure_pauli(ψ, N+1, "Z") # z-axis projection of ancilla qubit prepared in |+X> (google paper)

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
