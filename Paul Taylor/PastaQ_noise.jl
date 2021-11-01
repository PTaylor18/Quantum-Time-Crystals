using PastaQ
using ITensors
using Plots

# Define custom function to measure an observable, in this
# case a Pauli operator on `site`
function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

x_flip_circuit = [gatelayer("X", 1)]

σz1(ψ::MPS) = measure_pauli(ψ, 1, "Z")
σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]
σx(ψ::MPS) = [measure_pauli(ψ, j, "X") for j in 1:length(ψ)]

function noise_evolve(N, nsteps)
      t_vec = Vector{Float64}();
      Mz_vec = Vector{Float64}();
      circuit = Vector[];
      coupling_sequence = coupling_seq(N);

      layer = [gatelayer("X", 1)]

      obs = Observer([
            "χs" => linkdims,      # bond dimension at each bond
            "χmax" => maxlinkdim,  # maximum bond dimension
            "σᶻ(1)" => σz1,      # pauli Z on site 11
            "σᶻ" => σz,            # pauli Z on all sites
            ])

      for i=1:nsteps
          append!(t_vec, i)
          push!(circuit, layer)
      end

      noise = true
      if noise
          ψ = runcircuit(circuit; noise = ("amplitude_damping", (γ = 1,)), (observer!)=obs)
      else
          ψ = runcircuit(circuit; (observer!)=obs)
      end

      Mz_res = results(obs, "σᶻ(1)") # results for σᶻ
      for i=1:length(Mz_res)
          append!(Mz_vec, mean(Mz_res[i]))
      end

      return t_vec, Mz_vec

end

t_vec, Mz_vec = noise_evolve(20, 250);
plot(t_vec[1:100], Mz_vec[1:100], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",  legend=false)

N = 4
depth = 4
nshots = 5
layer = [("X", 1)]
circuit = randomcircuit(N, depth)

circuit = Vector[]

push!(circuit, layer)
circuit

# 1. Generation of measurement data on the quantum states
# at the output of a circuit. Each data-point is a projetive
# measurement in an arbitrary local basis. The default local basis
# is `["X","Y","Z"]`.
# a) Unitary circuit
# Returns output state as MPS
println("Generate samples from random projective measurements of the state U|0,0,…>:")
data, ψ = getsamples(
  circuit, nshots; local_basis=["X", "Y", "Z"], noise=("amplitude_damping", (γ=0.01,))
)
# Example of writing and reading
writesamples(data, ψ, "data/qst_circuit.h5")
data, ψ = readsamples("data/qst_circuit.h5")
@show maxlinkdim(ψ)
display(data)
println()
