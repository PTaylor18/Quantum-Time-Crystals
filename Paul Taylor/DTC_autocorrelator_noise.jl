using PastaQ
using ITensors
using Plots
using Statistics
using Distributions

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

function noise_correlator_evolve(N::Int, cor_pos::Int, nsteps::Int)
    t_vec = Vector{Float64}();
    cor_vec = Vector{Float64}();

    coupling_sequence = coupling_seq(N);
    circuit = Vector[];

    # first layer - can control which qubit is measured for auto-correlator measurement
    layer = [[("H", N+1)],[("CZ", (cor_pos, N+1))], gatelayer("Rx", N; (θ=π*97)),
            [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)],
            [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
            [("CZ", (cor_pos, N+1))],
            [("H", N+1)]]

    nshots = 20

    for i=1:nsteps
        append!(t_vec, i)
        push!(circuit, layer)

        data, ρ = getsamples(circuit, nshots; local_basis=["Z"], noise = ("amplitude_damping", (γ = 0.0001,)))
        #data, ρ = getsamples(circuit, nshots; local_basis=["Z"])

        res = Vector{Float64}();
        ancilla_res = data[:, cor_pos] # only look at ancilla qubit measurements

        for x in ancilla_res
            append!(res, 2 * (x.second - 0.5))
        end

        append!(cor_vec, mean(res))
        println("Step: $i")
    end
    return t_vec, cor_vec
end

nsteps = 100

t_vec, cor_vec = noise_correlator_evolve(20, 11, nsteps);
plot(t_vec[1:nsteps], cor_vec[1:nsteps], linetype=:steppre, xaxis="Time T", yaxis="<ψ|Ẑ(0)Ẑ(t)|ψ>",  legend=false)


# maybe do individual noise
# start with zero state then run circuit with noise with each layer
# try get no noise on ancilla qubit
# run circuit with cz no noise then run x str with noise zz with noise then cz no noise
# then use get samples




coupling_sequence = coupling_seq(N);
circuit = Vector[];
coupling_sequence_test = [(1,2),(2,3),(3,4)]

# first layer - can control which qubit is measured for auto-correlator measurement
layer = [[("H", 5)],[("CZ", (3, 5))], gatelayer("Rx", 4; (θ=π*0.97)),
        [("ZZ_couple", coupling_sequence_test[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence_test)],
        [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:4],
        [("CZ", (3, 5))],
        [("H", 5)]]

push!(circuit, layer)
circuit

data, ρ = getsamples(circuit, nshots; local_basis=["Z"])

ancilla_res = data[:,5]

res = Vector{Float64}();
for x in ancilla_res
    append!(res, 2 * (x.second - 0.5))
end

res
