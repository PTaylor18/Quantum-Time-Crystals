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

    layer = [gatelayer("Rx", 4; (θ=π*0.92)),
            [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)],
            [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]]

    nshots = 100

    for i=1:nsteps
        append!(t_vec, i)
        push!(circuit, layer)

        #data, ρ = getsamples(circuit, nshots; local_basis=["Z"], noise = ("dephasing", (γ = 0.1,)))
        data, ρ = getsamples(circuit, nshots; local_basis=["Z"])

        res = Vector{Float64}();

        for x in data
            append!(res, 2 * (x.second - 0.5))
        end

        append!(Mz_vec, mean(res))
        println("Step: $i")
    end
    return t_vec, Mz_vec
end

nsteps = 50

t_vec, Mz_vec = noise_evolve2(4, nsteps);
plot(t_vec[1:nsteps], Mz_vec[1:nsteps], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",  legend=false)


begin
    coupling_sequence = coupling_seq(N);

    layer = [gatelayer("Rx", 4; (θ=π*0.97)),
            [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)],
            [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]]

end

function getshots(circ, nshots, j)
    data, ψ = getsamples(circ, nshots; local_basis=["Z"], noise = ("amplitude_damping", (γ = 0.001,)))
    results = data[:, j]
    measurement_list = [];
    for n=1:nshots
        measured_bit = parse(Float64, replace("Z", results[n]))
        measured_spin = 2 * (measured_bit - 0.5)
        append!(measurement_list, measured_spin)
    end
    return measurement_list
end


nshots = 10
nsteps = 20
circuit = Vector[];

for i=1:nsteps
    push!(circuit, layer)
    measured_shots = getshots(circuit, nshots, 1)
    println(mean(measured_shots))
end

measured_shots
