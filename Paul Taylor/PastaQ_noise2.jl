using PastaQ
using ITensors
using Plots
using Statistics
using Distributions

obs = Observer([
      "χs" => linkdims,      # bond dimension at each bond
      "χmax" => maxlinkdim,  # maximum bond dimension
      "σᶻ(1)" => σz1,      # pauli Z on site 11
      "σᶻ" => σz,            # pauli Z on all sites
      ])

# pauli Z on all sites
σx2(ψ::MPS) = measure_pauli(ψ, 2, "X")
σz11(ψ::MPS) = measure_pauli(ψ, 11, "Z")
σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]


ρ = productstate(4)
gates = [gatelayer("X", 1)]

ρ = runcircuit(ρ, gates; (observer!)=obs)
ρ = runcircuit(ρ, gates; noise = ("amplitude_damping", (γ = 0.1,)))

nshots = 5
getsamples(ρ, nshots=5, local_basis = ["Z"])

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

function coupling_seq(N)
    sequence  = []
    for i=1:1:N-1
      push!(sequence, (i,i+1))
    end
    return sequence
end

function noise_evolve2(N, nsteps)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();

    coupling_sequence = coupling_seq(N);
    circuit = Vector[];

    layer = [gatelayer("Rx", 4; (θ=π*0.97)),
            [("ZZ_couple", coupling_sequence_test[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence_test)],
            [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:4]]

    nshots = 10

    for i=1:nsteps
        append!(t_vec, i)
        push!(circuit, layer)

        data = getsamples(circuit, nshots; local_basis=["Z"], noise = ("amplitude_damping", (γ = 0.00001,)))[1]

        res = Vector{Float64}();

        for x in data
            append!(res, x.second)
        end

        append!(Mz_vec, mean(res))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = noise_evolve2(10,250);
plot(t_vec[1:250], Mz_vec[1:250], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",  legend=false)
