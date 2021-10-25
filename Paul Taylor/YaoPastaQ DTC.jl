using YaoPastaQ, Yao
using PastaQ: MPS
using ITensors
using Random
using Plots
using Statistics
using Distributions

# use YaoPastaQ to convert a Yao circuit to a list of gates recognised by PastaQ

function yao_dtc_to_list(n)

    cirq = chain(n)
    push!(cirq, chain(n, prod([put(n, i=>Rx(π*0.97)) for i=1:n]))) # RxStr
    push!(cirq, prod([chain(n, put(n, i=>Rz(rand(Uniform(-π*1.5,-π*0.5))))*put(n, i+1=>Rz(rand(Uniform(-π*1.5,-π*0.5))))) for i=1:n-1])) # ZZpairs
    push!(cirq, chain(n, prod([put(n, i=>Rz(rand(Uniform(-π,π)))) for i=1:n]))) # RzStr

    list = []
    push!(list, genlist(chain(n, prod([put(n, i=>Rx(π*0.97)) for i=1:n]))))
    push!(list, genlist(prod([chain(n, put(n, i=>Rz(rand(Uniform(-π*1.5,-π*0.5))))*put(n, i+1=>Rz(rand(Uniform(-π*1.5,-π*0.5))))) for i=1:n-1])))
    push!(list, genlist(chain(n, prod([put(n, i=>Rz(rand(Uniform(-π,π)))) for i=1:n]))))

    return list
end


gates = [("X" , 1),                        # Pauli X on qubit 1
         ("CX", (1, 3)),                   # Controlled-X on qubits [1,3]
         ("Rx", 2, (θ=0.5,)),              # Rotation of θ around X
         ("Rn", 3, (θ=0.5, ϕ=0.2, λ=1.2)), # Arbitrary rotation with angles (θ,ϕ,λ)
         ("√SWAP", (3, 4)),                # Sqrt Swap on qubits [2,3]
         ("T" , 4)]

cirq = yao_dtc_to_list(4)
cirq = test_circuit = [gatelayer("Rx", 4; (θ=π*0.97))]
obs = Observer([
  "χs" => linkdims,      # bond dimension at each bond
  "χmax" => maxlinkdim,  # maximum bond dimension
  "σᶻ" => σz,            # pauli Z on all sites
  ])


runcircuit(gates; (observer!)=obs)

function Mz_evolve(N, nsteps)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    Mz_11_vec = Vector{Float64}();
    circuit = Vector[];
    # first layer
    layer = yao_dtc_to_list(N)
    for i=1:nsteps
        append!(t_vec, i)
        push!(circuit, layer)
    end

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

    noise = false
    if noise
        ψ = runcircuit(circuit; (observer!)=obs, noise = ("amplitude_damping", (γ = 0.01,)))
    else
        ψ = runcircuit(circuit; (observer!)=obs)
    end

    ψ = runcircuit(circuit; (observer!)=obs)

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

t_vec, Mz_vec, Mz_11_vec = Mz_evolve(20, 250);
plot(t_vec[1:250], Mz_vec[1:250], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",  legend=false)
