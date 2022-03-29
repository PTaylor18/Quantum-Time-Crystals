using PastaQ
using ITensors
using Random
using Plots
using Statistics
using Distributions
using FFTW
using JLD
using DataFrames
using CSV
using DelimitedFiles
using StatsBase
using LaTeXStrings
using LinearAlgebra

import PastaQ: gate

# Ising (ZZ) coupling gate
function gate(::GateName"ZZ"; ϕ::Number)
  return [
    exp(-im*ϕ/4) 0 0 0
    0 exp(im*ϕ/4) 0 0
    0 0 exp(im*ϕ/4) 0
    0 0 0 exp(-im*ϕ/4)
  ]
end

gate(::GateName"ZZ_couple"; kwargs...) = gate("ZZ"; kwargs...)

function disorder_arrays(N, r, nsteps)

    # prepare vector of arrays for multiple disorder instances
    j_run_disorder_array = Array[]
    h_run_disorder_array = Array[]

    for k=1:r

        j_disorder_array = zeros(1, N + (-1)^(4%3))
        h_disorder_array = zeros(1, N)

        for i=1:nsteps

            randϕ = rand(Uniform(-0.375π, -0.125π), N + (-1)^(4%3)) # j disorder
            j_disorder_array = vcat(j_disorder_array, randϕ')
        end

        for i=1:nsteps
            randθ = rand(Uniform(-π,π), N) # longitudinal field disorder
            h_disorder_array = vcat(h_disorder_array, randθ')
        end

        j_disorder_array = j_disorder_array[Not(1),:]
        h_disorder_array = h_disorder_array[Not(1),:]

        push!(j_run_disorder_array, j_disorder_array)
        push!(h_run_disorder_array, h_disorder_array)
    end
    return j_run_disorder_array, h_run_disorder_array
end

disorder_instance = disorder_arrays(N, r, 60)

j_disorder_array = disorder_instance[1][1]
h_disorder_array = disorder_instance[2][k]

r =5
test = 1

#---

function spin_glass(N, r, initial_state::String="polarized")

  g_vec = collect(0.7:0.02:1.0)
  χ_sg_final = Vector{Float64}(); # final average spin glass for each g

  bond_dim_final = Vector{Float64}(); # final average bond dim for each g

  coupling_sequence = PastaQ.lineararray(N)

  # generate an initial state
  function initial_state_gen(N, initial_state)

      if initial_state == "polarized"   # initalise in |000...> state

          pol_vec = convert(Vector{Int64}, zeros(N))
          bit_string = prod(string.(pol_vec))
          pol_ancilla = append!(pol_vec, 0)
          pol_state = [("I", i) for i in 1:N+1]
          ini_state = pol_state

      elseif initial_state == "neel"    # initialise in |01010...> state

          n = convert(Vector{Int64}, ones(N))
          neel_string = [isodd(i) == true ? 0 : 1 for i in 1:N]
          bit_string = prod(string.(neel_string))
          neel_ancilla = append!(neel_string, 0)
          neel_state = [neel_ancilla[j] == 1 ? ("X", j) : ("I", j) for j in 1:N+1]
          ini_state = neel_state

      elseif initial_state == "random"     # initalise in |0111010010...> state

          random_bitstring = [parse(Int, x) for x in rand(("0", "1"), N)] # generates a random bistring
          bit_string = prod(string.(random_bitstring))
          bitstring_ancilla = append!(random_bitstring, 0) # add 0 for ancilla qubit
          random_bitstring_state = [bitstring_ancilla[n] == 1 ? ("X", n) : ("I", n) for n in 1:N+1] # generates bitstring circuit
          ini_state = random_bitstring_state
      end

      return ini_state, bit_string
  end

  ini_state, bit_string = initial_state_gen(N, initial_state)

  disorder_instance = disorder_arrays(N, r, 60)

  for g=0.7:0.02:1.0

      av_χ_sg = Vector{Float64}(); # average spin glass values for each run

      run_vec = Vector{Float64}(); # vector of mean spin glass for each run

      max_bond_dim_vec = Vector{Float64}();

      for k=1:r

          println("\nRunning g = $g with Disorder Instance $k")
          ψ = runcircuit(ini_state) # initialise MPS with initial state
          χ_sg_vec = Vector{Float64}(); # Vector of spin glass values between t=51 and t=60

          j_disorder_array = disorder_instance[1][k]
          h_disorder_array = disorder_instance[2][k]

          av_bond_dim = Vector{Float64}(); # average bond dim values for each run

          for i=1:60

              println("Step $i")

              # select disorder vectors for j and h for current step
              randϕ = j_disorder_array[i,:]
              randθ = h_disorder_array[i,:]

              # q controls which qubit is measured for auto-correlator measurement

              layer = [gatelayer("Rx", N; (θ=π*g)),
                      [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
                      [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                      [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                      ]

              # layer = [ gatelayer("Rx", N; (θ=π*g)),
              #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=-0.4,)) for i=1:length(coupling_sequence[1])],
              #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=-0.4,)) for i=1:length(coupling_sequence[2])],
              #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
              #         ]

              # layer = [gatelayer("Rx", N; (θ=π*g)),
              #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=0,)) for i=1:length(coupling_sequence[1])],
              #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=0,)) for i=1:length(coupling_sequence[2])],
              #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
              #         ]

              ψ = runcircuit(ψ, layer)

              # calculates χˢᵍ order parameter as in google paper

              zzcorr = correlation_matrix(ψ, "Z", "Z"; site_range=2:N-1)
              zzcorr_squared = zzcorr.^2
              χ_sg = 1/(N-2)*sum(tril(zzcorr_squared, -1)) # sum the lower triangle of the matrix

              append!(χ_sg_vec, real(χ_sg))

              append!(av_bond_dim, maxlinkdim(ψ))

          end

          χ_sg_mean = mean(χ_sg_vec[51:60]) # average between t=51 and 5=60
          println(χ_sg_mean)
          append!(av_χ_sg, χ_sg_mean)

          append!(max_bond_dim_vec, maximum(av_bond_dim))

      end
      append!(χ_sg_final, mean(trim(av_χ_sg, prop=0.3)))

      append!(bond_dim_final, mean(max_bond_dim_vec))
  end
  return g_vec, χ_sg_final, bond_dim_final
end

#---
# Simulations

r = 10
e
g_vec, χ_sg_final_8 = spin_glass(8, 10, "polarized");
g_vec, χ_sg_final_10 = spin_glass(10, r, "polarized");
g_vec, χ_sg_final_12 = spin_glass(12, r, "polarized");
g_vec, χ_sg_final_14 = spin_glass(14, r, "polarized");

g_vec, χ_sg_final, bond_dim_final = spin_glass(8, 5, "polarized");

plot(g_vec[1:16], bond_dim_final[1:16], xaxis=L"\textbf{g}", yaxis=L"\textbf{Spin \ Glass \ Order \ Parameter, \chi^{sg}}",
    label="N=8",
    marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:bottomright)


plot(g_vec[1:11], χ_sg_final_8[1:11], xaxis=L"\textbf{g}", yaxis=L"\textbf{Spin \ Glass \ Order \ Parameter, \chi^{sg}}",
    label="N=8",
    marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:topleft)
plot!(g_vec[1:11], χ_sg_final_10[1:11], marker=(:circle), label="N=10")
plot!(g_vec[1:11], χ_sg_final_12[1:11], marker=(:circle), label="N=12")
plot!(g_vec[1:11], χ_sg_final_14[1:11], marker=(:circle), label="N=14")

# saving data
begin
    save("sg_8_new.jld", "g_vec", g_vec, "sg_8_new", χ_sg_final_8)
    save("sg_10_new.jld", "g_vec", g_vec, "sg_10_new", χ_sg_final_10)
    save("sg_12_new.jld", "g_vec", g_vec, "sg_12_new", χ_sg_final_12)
    save("sg_14_new.jld", "g_vec", g_vec, "sg_14_new", χ_sg_final_14)
end

function data_to_array()

    χ_sg_final_8 = get(load("spinglass_data/sg_8.jld"), "sg_8", 1)
    χ_sg_final_10 = get(load("spinglass_data/sg_10.jld"), "sg_10", 1)
    χ_sg_final_12 = get(load("spinglass_data/sg_12.jld"), "sg_12", 1)
    χ_sg_final_14 = get(load("spinglass_data/sg_14.jld"), "sg_14", 1)

    A = [χ_sg_final_8,
        χ_sg_final_10,
        χ_sg_final_12,
        χ_sg_final_14]

    return A
end
g_vec[1:11]
dataarray = [χ_sg_final_8[1:11],
    χ_sg_final_10[1:11],
    χ_sg_final_12[1:11],
    χ_sg_final_14[1:11]]


writedlm("sg_new_data.csv", dataarray, ',')

#manual fss
β_sg = 1.07
ν_sg = 0.83
gc_sg = 0.88

# plot the scaled curves
begin
    # plot((g_vec[1:11].-gc_sg).*8^(1/ν_sg), χ_sg_final_8[1:11].*8^(-β_sg/ν_sg),
    # xlabel=L"(g-g_c)*N^{1/\nu}", ylabel=L"\chi^{sg} * N^{-\beta/\nu}",
    # label="N=8",
    # line=(:dash), marker=(:circle), palette= :seaborn_bright,
    # grid=false, legend=:topleft)
    plot((g_vec[1:11].-gc_sg).*10^(1/ν_sg), χ_sg_final_10[1:11].*10^(-β_sg/ν_sg),
    xlabel=L"(g-g_c)*N^{1/\nu}", ylabel=L"\chi^{sg} * N^{-\beta/\nu}",
    label="N=10",
    line=(:dash), marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:topleft)
    plot!((g_vec[1:11].-gc_sg).*12^(1/ν_sg), χ_sg_final_12[1:11].*12^(-β_sg/ν_sg), label="N=12", line=(:dash), marker=(:circle))
    plot!((g_vec[1:11].-gc_sg).*14^(1/ν_sg), χ_sg_final_14[1:11].*14^(-β_sg/ν_sg), label="N=14", line=(:dash), marker=(:circle))
    #plot!((g_vec[1:11].-gc_sg).*14^(1/ν_sg), χ_sg_final_14[1:11].*14^(-β_sg/ν_sg), label="N=14", line=(:dash),marker=(:circle))
    annotate!([(-4.5, 0.1, (latexstring("\\beta = $(β_sg)"), 14)),(-4.49, 0.082, (latexstring("\\nu = $(ν_sg)"), 14)), (-4.55, 0.060,(latexstring("g_c = $(gc_sg)"), 14))])
end


g_vec, χ_sg_final_8_nodis = spin_glass(8, r, "polarized");
g_vec, χ_sg_final_10_nodis = spin_glass(10, r, "polarized");
g_vec, χ_sg_final_12_nodis = spin_glass(12, r, "polarized");
g_vec, χ_sg_final_14_nodis = spin_glass(14, r, "polarized");

plot(g_vec[1:16], χ_sg_final_8_nodis[1:16], xaxis=L"\textbf{g}", yaxis=L"\textbf{Spin \ Glass \ Order \ Parameter, \chi^{sg}}",
    label="N=8",
    marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:topleft)
plot!(g_vec[1:16], χ_sg_final_10_nodis[1:16], marker=(:circle), label="N=10")
plot!(g_vec[1:16], χ_sg_final_12_nodis[1:16], marker=(:circle), label="N=12")

#---

# [-1.125π, -0.375π]
g_vec, χ_sg_final_8_34 = spin_glass(8, 20, "polarized");
g_vec, χ_sg_final_10_34 = spin_glass(10, r, "polarized");
g_vec, χ_sg_final_12_34 = spin_glass(12, r, "polarized");

plot(g_vec[1:12], χ_sg_final_8_34[1:12], xaxis=L"\textbf{g}", yaxis=L"\textbf{Spin \ Glass \ Order \ Parameter, \chi^{sg}}",
    label="N=8",
    marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:topleft)
plot!(g_vec[1:12], χ_sg_final_10_34[1:12], marker=(:circle), label="N=10")
plot!(g_vec[1:12], χ_sg_final_12_34[1:12], marker=(:circle), label="N=12")


# [-0.75π, -0.25π]
g_vec, χ_sg_final_8_12 = spin_glass(8, 20, "polarized");
g_vec, χ_sg_final_10_12 = spin_glass(10, r, "polarized");
g_vec, χ_sg_final_12_12 = spin_glass(12, r, "polarized");

plot(g_vec[1:13], χ_sg_final_8_12[1:13], xaxis=L"\textbf{g}", yaxis=L"\textbf{Spin \ Glass \ Order \ Parameter, \chi^{sg}}",
    label="N=8",
    marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:topleft)
plot!(g_vec[1:13], χ_sg_final_10_12[1:13], marker=(:circle), label="N=10")
plot!(g_vec[1:13], χ_sg_final_12_12[1:13], marker=(:circle), label="N=12")


# [-0.375π, -0.125π]
g_vec, χ_sg_final_8_14 = spin_glass(8, 20, "polarized");
g_vec, χ_sg_final_10_14 = spin_glass(10, r, "polarized");
g_vec, χ_sg_final_12_14 = spin_glass(12, r, "polarized");

plot(g_vec[1:13], χ_sg_final_8_14[1:13], xaxis=L"\textbf{g}", yaxis=L"\textbf{Spin \ Glass \ Order \ Parameter, \chi^{sg}}",
    label="N=8",
    marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:topleft)
plot!(g_vec[1:13], χ_sg_final_10_14[1:13], marker=(:circle), label="N=10")
plot!(g_vec[1:13], χ_sg_final_12_14[1:13], marker=(:circle), label="N=12")


# with no interactions
g_vec, χ_sg_final_8_0 = spin_glass(8, 10, "polarized");
g_vec, χ_sg_final_10_0 = spin_glass(10, r, "polarized");
g_vec, χ_sg_final_12_0 = spin_glass(12, r, "polarized");

plot(g_vec[8:16], χ_sg_final_8_0[8:16], xaxis=L"\textbf{g}", yaxis=L"\textbf{Spin \ Glass \ Order \ Parameter, \chi^{sg}}",
    label="N=8",
    marker=(:circle), palette= :seaborn_bright,
    grid=false, legend=:topleft)
plot!(g_vec[8:16], χ_sg_final_10_0[8:16], marker=(:circle), label="N=10")
plot!(g_vec[1:11], χ_sg_final_12_0[1:11], marker=(:circle), label="N=12")


N=8


coupling_sequence = PastaQ.lineararray(N)

ψ = productstate(N)

bonddimvec = Vector{Float64}()
for i=1:30
    layer = [gatelayer("Rx", 4; (θ=π*0.7)),
            [("ZZ_couple", coupling_sequence[1][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[1])],
            [("ZZ_couple", coupling_sequence[2][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[2])],
            [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]]
    ψ = runcircuit(ψ, layer)
    append!(bonddimvec, maxlinkdim(ψ))
end

bonddimvec

linkdim(ψ, 3)
ψ
maxlinkdim(ψ)

PastaQ.measure()

ψ = runcircuit(ψ, layer)

measure(ψ, (linkdim) )
