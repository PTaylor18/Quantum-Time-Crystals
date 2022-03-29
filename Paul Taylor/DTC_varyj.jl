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

Z_field(N::Int) = [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]

# pauli Z on all sites

σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]



#---

function Mz_varyj_evolve(N, g, j, σ, ϵ, nsteps)

    """
    evolves a DTC with an average j disorder and variance σ

    """

    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();

    χ_sg_final = Vector{Float64}(); # final average spin glass for each g

    coupling_sequence = PastaQ.lineararray(N)
    circuit = Vector[];

    neel_state = [("X", i) for i=2:2:N]
    #push!(circuit, neel_state)
    #circuit = randombitstring(N)

    #circuit = push!(circuit, [("X", 1),("X", 3),("X", 4),("X", 7)])

    ψ = productstate(N)
    #ψ = runcircuit(ψ, neel_state)

    # define the circuit observer

    obs = Observer([
      "χs" => linkdims,      # bond dimension at each bond
      "χmax" => maxlinkdim,  # maximum bond dimension
      "σᶻ" => σz,            # pauli Z on all sites
      ])


    for i=1:nsteps
        append!(t_vec, i)

        if j == 0 # no interactions
            layer = [gatelayer("Rx", N; (θ=(π*g))),
                    [("ZZ_couple", coupling_sequence[1][i], (ϕ=0,)) for i=1:length(coupling_sequence[1])],
                    [("ZZ_couple", coupling_sequence[2][i], (ϕ=0,)) for i=1:length(coupling_sequence[2])],
                    [("Rz", i, (ϕ=rand(Uniform(0,π)),)) for i=1:N]]
        else
            randϕ = rand(Normal(j, σ), N + (-1)^(4%3)) # set disordered angles for Ising coupling

            # layer = [gatelayer("Rx", N; (θ=π*g)),
            #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
            #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
            #         [("Rz", i, (ϕ=rand(Uniform(0,π)),)) for i=1:N]]

            layer = [gatelayer("Rx", N; (θ=(π*g))),
                    [("ZZ_couple", coupling_sequence[1][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[1])],
                    [("ZZ_couple", coupling_sequence[2][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[2])],
                    [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]]
            # rand(Uniform(-π,π))
        end

        randθ = rand(Normal(π, σ), N)

        push!(circuit, layer)
    end


    ψ = runcircuit(ψ, circuit; (observer!)=obs)

    # obtain magnetizations for all qubits at all evolutions
    Mz_res = results(obs, "σᶻ") # results for σᶻ
    for i=1:length(Mz_res)
        append!(Mz_vec, mean(Mz_res[i]))
    end

    # create array of magnetizations for each qubit
    mz_q_array = zeros(1,nsteps)

    for q=1:N
        mz_q_vec = Vector{Float64}()
        for i=1:nsteps
            append!(mz_q_vec, Mz_res[i][q])
        end
        mz_q_array = vcat(mz_q_array, mz_q_vec')
    end

    mz_q_array = mz_q_array[Not(1),:]

    return t_vec, Mz_vec, mz_q_array
end

N = 14 # No. of qubits
g = 0.85
g = π
j = 0.25
j = -1.571
σ = abs.(0.2*j)
σ = 0.1
nsteps = 100
ϵ = 0.2


π - ϵ
randϕ = rand(Normal(j, σ), N + (-1)^(4%3))
randϕ = zeros(N + (-1)^(4%3))

coupling_sequence = PastaQ.lineararray(N)

t_vec, Mz_vec, mzt_q_array = Mz_varyj_evolve(N, g, j, σ, ϵ, nsteps);

# plot average magnetization
begin
    plot(t_vec[1:nsteps], Mz_vec[1:nsteps],
    linetype=:steppre, xaxis=L"t", yaxis=L"Average \ Magnetisation \ Mz", ylims=(-1.05,1),
    label=latexstring("N=$N, g=$g"), title=L"Thermal", grid=false, legend=:topright)
end
# can plot magnetization of each individual spin
begin
    q = 11 # select qubit
    plot(t_vec[1:nsteps], mzt_q_array[q,:][1:nsteps],
        linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",
        label="Qubit $q, g=$g", legend=:topright)
end

Mz_vec[10:nsteps]
# plot fourier transform
begin
    Mz_vec_fft = fft(Mz_vec[10:nsteps])
    freq_vec = (t_vec[10:nsteps].-10)./(nsteps-10)
    plot(freq_vec, 1/maximum(abs.(Mz_vec_fft))*abs.(Mz_vec_fft), xaxis=L"\omega/\omega_D", yaxis =L"FFT", grid=false, label=label=latexstring("g=$(g)"), legend=:topright)
    plot!(size=(300,200))
end


L"(g-g_c)*N^{1/\nu}"

plot()

inset_sub

d = Normal(j, abs.(0.2*j))
plot(x->pdf(d,x), xaxis = "x", yaxis="Probability Density", legend=false)

function randombitstring(N)

    circuit = Vector[]

    randstring = vec(rand(1:N, N,1))

    bitstring = [("X", i) for i in randstring]

    push!(circuit, bitstring)

end

#---

# testing without loop

Mz_vec = Vector{Float64}();
coupling_sequence = PastaQ.lineararray(N)
circuit = Vector[];

neel_state = [("X", i) for i=2:2:N]
ψ = productstate(N)
ψ = runcircuit(ψ, neel_state)

obs = Observer([
  "χs" => linkdims,      # bond dimension at each bond
  "χmax" => maxlinkdim,  # maximum bond dimension
  "σᶻ" => σz,            # pauli Z on all sites
  ])
for i=1:nsteps
    layer = [gatelayer("Rx", N; (θ=π*g)),
            [("ZZ_couple", coupling_sequence[1][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[1])],
            [("ZZ_couple", coupling_sequence[2][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[2])],
            [("Rz", i, (ϕ=rand(Uniform(0,π)),)) for i=1:N]]

    push!(circuit, layer)

end

ψ = runcircuit(ψ, circuit; (observer!)=obs)

Mz_res = results(obs, "σᶻ") # results for σᶻ
for i=1:length(Mz_res)
    append!(Mz_vec, mean(Mz_res[i]))
end

Mz_vec
plot(t_vec[1:nsteps], Mz_vec[1:nsteps],
    linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",
    label="N=$N, g=$g, J=$j, σ=$σ", legend=:topright)
