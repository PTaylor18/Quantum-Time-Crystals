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
    #push!(cirq, prod([chain(n, put(n, i=>Rz(rand(Uniform(-π*1.5,-π*0.5))))*put(n, i+1=>Rz(rand(Uniform(-π*1.5,-π*0.5))))) for i=1:n-1])) # ZZpairs
    push!(cirq, prod([chain(n, put(n, i=>Z)*put(n, i+1=>Z)) for i=1:n-1])) # ZZpairs
    push!(cirq, chain(n, prod([put(n, i=>Rz(rand(Uniform(-π,π)))) for i=1:n]))) # RzStr with disorder

    list = Vector[]
    push!(list, genlist(chain(n, prod([put(n, i=>Rx(π*1.)) for i=1:n]))))
    push!(list, genlist(prod([chain(n, put(n, i=>Z)*put(n, i+1=>Z)) for i=1:n-1])))
    push!(list, genlist(chain(n, prod([put(n, i=>Rz(rand(Uniform(-π,π)))) for i=1:n]))))

    return list
end

yao_dtc_to_list(4)

function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

σx2(ψ::MPS) = measure_pauli(ψ, 2, "X")
σz11(ψ::MPS) = measure_pauli(ψ, 11, "Z")
σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]

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
plot(t_vec[1:50], Mz_vec[1:50], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",  legend=false)
