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

function gate(::GateName"XX"; ϕ::Number)
  return [
    cos(ϕ/2) 0 0 -im*sin(ϕ/2)
    0 cos(ϕ/2) -im*sin(ϕ/2) 0
    0 -im*sin(ϕ/2) cos(ϕ/2) 0
    -im*sin(ϕ/2) 0 0 cos(ϕ/2)
  ]
end

gate(::GateName"XX_couple"; kwargs...) = gate("XX"; kwargs...)

function gate(::GateName"YY"; ϕ::Number)
  return [
    cos(ϕ/2) 0 0 im*sin(ϕ/2)
    0 cos(ϕ/2) -im*sin(ϕ/2) 0
    0 -im*sin(ϕ/2) cos(ϕ/2) 0
    im*sin(ϕ/2) 0 0 cos(ϕ/2)
  ]
end

gate(::GateName"YY_couple"; kwargs...) = gate("YY"; kwargs...)

function gate(::GateName"XY"; ϕ::Number)

    if ϕ == -π || π
        return [
          1 0 0 0
          0 0 im 0
          0 im 0 0
          0 0 0 1
        ]
    else
        return [
          1 0 0 0
          0 cos(ϕ/2) -im*sin(ϕ/2) 0
          0 -im*sin(ϕ) cos(ϕ/2) 0
          0 0 0 1
        ]
    end
end

function gate(::GateName"XY"; ϕ::Number)
    return [
      1 0 0 0
      0 cos(ϕ/2) -im*sin(ϕ/2) 0
      0 -im*sin(ϕ) cos(ϕ/2) 0
      0 0 0 1
    ]
end

gate(::GateName"XY_couple"; kwargs...) = gate("XY"; kwargs...)

function gate(::GateName"XX+YY"; ϕ::Number)

  if ϕ == -π/2
    return [
      1 0 0 0
      0 0 im 0
      0 im 0 0
      0 0 0 1
    ]
  else
  return
end

gate(::GateName"XX+YY_couple"; kwargs...) = gate("XX+YY"; kwargs...)



a == b
rand(Uniform(-π/2,-π/2+0.1))
function disorder_arrays(N, r, nsteps)

    # prepare vector of arrays for multiple disorder instances
    j_run_disorder_array = Array[]
    h_run_disorder_array = Array[]
    xy_run_disorder_array = Array[]

    for k=1:r

        j_disorder_array = zeros(1, N + (-1)^(4%3))
        xy_disorder_array = zeros(1, N + (-1)^(4%3))
        h_disorder_array = zeros(1, N)

        for i=1:nsteps

            randϕ = rand(Uniform(-π*1.5,-π*0.5), N + (-1)^(4%3)) # j disorder
            j_disorder_array = vcat(j_disorder_array, randϕ')

            randω = rand(Uniform(-π*1.5,-π*0.5), N + (-1)^(4%3)) # xy disorder
            xy_disorder_array = vcat(xy_disorder_array, randω')
        end

        for i=1:nsteps
            randθ = rand(Uniform(-π,π), N) # longitudinal field disorder
            h_disorder_array = vcat(h_disorder_array, randθ')
        end

        j_disorder_array = j_disorder_array[Not(1),:]
        xy_disorder_array = xy_disorder_array[Not(1),:]
        h_disorder_array = h_disorder_array[Not(1),:]

        push!(j_run_disorder_array, j_disorder_array)
        push!(xy_run_disorder_array, xy_disorder_array)
        push!(h_run_disorder_array, h_disorder_array)
    end
    return j_run_disorder_array, xy_run_disorder_array, h_run_disorder_array
end

disorder_instance = disorder_arrays(N, r, 60)

j_disorder_array = disorder_instance[1][1]
xy_disorder_array = disorder_instance[2][1]
h_disorder_array = disorder_instance[3][1]

r =5
test = 1

#---

function XY_evolve(N::Int, g::Float64, j::Float64, σ::Float64,
    nsteps::Int, r::Int, disorder_instance=false, initial_state::String="polarized")

    """
    evolves a DTC with an average j disorder and variance σ and calculates the correlation function:

        <ψ|Ẑ(0)Ẑ(t)|ψ>

        N = no. of qubits
        g = X flip perturbation
        nsteps = number of discrete time evolutions
        r = number of runs to average on
        constant_disorder = if true then use a predefined disorder instance, if false then generate new array of disorders
        initial_state: select between polarized |0>^⊗N, neel |01>^⊗N/2 or random |0111010010...> states

    """

    t_vec= collect(1:nsteps)

    av_correlator_array = zeros(1,nsteps)

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

    final_av_correlator_vec = Vector{Float64}();

    for k=1:r

        println("\nRunning Disorder Instance $k")

        # generate disorder instances
        if disorder_instance == true
            println("Using predefined disorder instance")
            j_disorder_array = predefined_disorder_instance[1][k]
            xy_disorder_array = predefined_disorder_instance[2][k]
            h_disorder_array = predefined_disorder_instance[3][k]
        elseif disorder_instance == false
            full_array = disorder_arrays(N, 1, nsteps)
            j_disorder_array = full_array[1][1]
            h_disorder_array = full_array[2][1]
        end

        correlator_array = zeros(1,nsteps)
        av_auto_correlator_vec = Vector{Float64}();

        for q=2:1:N-1 # ignore edge qubits

            Mz_auto_correlator_res = Vector{Float64}();

            ψ = runcircuit(ini_state) # initialise MPS with initial state

            println("Running autocorrelator on qubit $q")

            for i=1:nsteps

                # select disorder vectors for j and h for current step
                randϕ = j_disorder_array[i,:]
                randω = xy_disorder_array[i,:]
                randθ = h_disorder_array[i,:]

                # q controls which qubit is measured for auto-correlator measurement

                # layer = [[("H", N+1)],[("CZ", (N+1, q))],
                #         gatelayer("Rx", N; (θ=π*g)),
                #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
                #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                #         [("CZ", (N+1, q))],
                #         [("H", N+1)]]

                # layer = [[("H", N+1)],[("CZ", (N+1, q))],
                #         gatelayer("Rx", N; (θ=π*g)),
                #         [("iSwap", coupling_sequence[1][i]) for i=1:length(coupling_sequence[1])],
                #         [("iSwap", coupling_sequence[2][i]) for i=1:length(coupling_sequence[2])],
                #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                #         [("CZ", (N+1, q))],
                #         [("H", N+1)]]

                layer = [[("H", N+1)],[("CZ", (N+1, q))],
                        gatelayer("Rx", N; (θ=π*g)),
                        [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
                        [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                        [("iSwap", coupling_sequence[1][i]) for i=1:length(coupling_sequence[1])],
                        [("iSwap", coupling_sequence[2][i]) for i=1:length(coupling_sequence[2])],
                        [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                        [("CZ", (N+1, q))],
                        [("H", N+1)]]

                # layer = [[("H", N+1)],[("CZ", (N+1, q))],
                #         gatelayer("Rx", N; (θ=π*g)),
                #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
                #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                #         [("XX+YY_couple", coupling_sequence[1][i], (ϕ=randω[i],)) for i=1:length(coupling_sequence[1])],
                #         [("XX+YY_couple", coupling_sequence[2][i], (ϕ=randω[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                #         [("CZ", (N+1, q))],
                #         [("H", N+1)]]

                # layer = [[("H", N+1)],[("CZ", (N+1, q))],
                #         gatelayer("Rx", N; (θ=π*g)),
                #         [("YY_couple", coupling_sequence[1][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[1])],
                #         [("YY_couple", coupling_sequence[2][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[2])],
                #         [("XX_couple", coupling_sequence[1][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[1])],
                #         [("XX_couple", coupling_sequence[2][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[2])],
                #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                #         [("CZ", (N+1, q))],
                #         [("H", N+1)]]

                # layer = [[("H", N+1)],[("CZ", (N+1, q))],
                #         gatelayer("Rx", N; (θ=π*g)),
                #         [("XY_couple", coupling_sequence[1][i], (ϕ=randω[i],)) for i=1:length(coupling_sequence[1])],
                #         [("XY_couple", coupling_sequence[2][i], (ϕ=randω[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
                #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                #         [("CZ", (N+1, q))],
                #         [("H", N+1)]]

                layer = [[("H", N+1)],[("CZ", (N+1, q))],
                        gatelayer("Rx", N; (θ=π*g)),
                        [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
                        [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                        [("XY_couple", coupling_sequence[1][i], (ϕ=randω[i],)) for i=1:length(coupling_sequence[1])],
                        [("XY_couple", coupling_sequence[2][i], (ϕ=randω[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                        [("CZ", (N+1, q))],
                        [("H", N+1)]]

                # layer = [[("H", N+1)],[("CZ", (N+1, q))],
                #         gatelayer("Rx", N; (θ=π*g)),
                #         [("XX+YY_couple", coupling_sequence[1][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[1])],
                #         [("XX+YY_couple", coupling_sequence[2][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[2])],
                #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                #         [("CZ", (N+1, q))],
                #         [("H", N+1)]]

                # layer = [[("H", N+1)],[("CZ", (N+1, q))],
                #         gatelayer("Rx", N; (θ=π*g)),
                #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=-0.4,)) for i=1:length(coupling_sequence[1])],
                #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=-0.4,)) for i=1:length(coupling_sequence[2])],
                #         [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                #         [("CZ", (N+1, q))],
                #         [("H", N+1)]]

                ψ = runcircuit(ψ, layer)

                Z = PastaQ.measure(ψ, ("Z", N+1)) # measures <Z> on the ancilla qubit
                append!(Mz_auto_correlator_res, Z)
            end

            # results for <ψ|Ẑ(0)Ẑ(t)|ψ> on control qubit q

            correlator_array = vcat(correlator_array, Mz_auto_correlator_res')
        end

        correlator_array = correlator_array[Not(1),:]
        av_auto_correlator_vec = vec(mean(correlator_array, dims=1))

        av_correlator_array = vcat(av_correlator_array, av_auto_correlator_vec')

    end

    av_correlator_array = av_correlator_array[Not(1),:]
    final_av_correlator_vec = vec(mean(av_correlator_array, dims=1))

    println("\nInitial bitstring: ", bit_string)

    return t_vec, av_correlator_array, final_av_correlator_vec, bit_string
end

N = 8 # No. of qubits
g = 0.94
r = 5
j = -0.15
σ = abs.(0.2*j)
σ = 0.2
nsteps = 50

initial_state = "random";

t_vec = collect(1:100)

predefined_disorder_instance = disorder_arrays(N, r, nsteps);
t_vec, correlator_array_pol, final_av_correlator_vec_pol, bit_string_pol = XY_evolve(N, g, j, σ, nsteps, r, true, "polarized");
t_vec, correlator_array_neel, final_av_correlator_vec_neel, bit_string_neel = XY_evolve(N, g, j, σ, nsteps, r, true, "neel");
t_vec, correlator_array_rand, final_av_correlator_vec_rand, bit_string_rand = XY_evolve(N, g, j, σ, nsteps, r, true, "random");

zz_vec = final_av_correlator_vec_pol
iswap_vec = final_av_correlator_vec_pol
zzxy_vec = final_av_correlator_vec_pol

plot(t_vec[1:nsteps], zz_vec[1:nsteps],
    xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
    label="ZZ", title=latexstring("N=$N, g=$g"),
    marker=(:circle), markersize=5, color=:blueviolet, linewidth=1, grid=false, legend=:topright)
plot!(t_vec[1:nsteps], iswap_vec[1:nsteps],
    xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
    label="ZZ + iSwap", title=latexstring("N=$N, g=$g"),
    marker=(:circle), markersize=5, color=:red3, linewidth=1, grid=false, legend=:topright)
plot!(t_vec[1:nsteps], zzxy_vec[1:nsteps],
    xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
    label="ZZ + XY", title=latexstring("N=$N, g=$g"),
    marker=(:circle), markersize=5, color=:black, linewidth=1, grid=false, legend=:topright)

plot(t_vec[1:nsteps], final_av_correlator_vec_pol[1:nsteps],
    xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
    label="ZZ", title=latexstring("N=$N, g=$g"),
    marker=(:circle), markersize=5, color=:blueviolet, linewidth=1, grid=false, legend=:topright)

begin
    plot(t_vec[1:nsteps], final_av_correlator_vec_pol[1:nsteps],
        xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
        label="Polarized", title=latexstring("N=$N, g=$g"),
        marker=(:circle), markersize=5, color=:blueviolet, linewidth=1, grid=false, legend=:topright)
    plot!(t_vec[1:nsteps], final_av_correlator_vec_neel[1:nsteps],
        label="Néel", marker=(:circle), markersize=5, linewidth=1.5, color=:red3)

    plot!(t_vec[1:nsteps], final_av_correlator_vec_rand[1:nsteps],
        label="Random", marker=(:circle), markersize=5, linewidth=2.5, color=:black)
end


begin
    ac_pol_vec_fft = normalize!(fft(final_av_correlator_vec_pol[10:nsteps]))
    freq_vec = (t_vec[10:nsteps].-10)./(nsteps-10)
    plot(freq_vec, (abs.(ac_pol_vec_fft)), xaxis=L"\omega/\omega_D", yaxis=L"FFT",
        label=latexstring("g = 1.0"), color=:red3, grid=false, legend=:topright)
end

layer = [[("H", N+1)],[("CZ", (N+1, 4))],
        gatelayer("Rx", N; (θ=π*g)),
        [("iSwap", coupling_sequence[1][i]) for i=1:length(coupling_sequence[1])],
        [("iSwap", coupling_sequence[2][i]) for i=1:length(coupling_sequence[2])],
        [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
        [("CZ", (N+1, 4))],
        [("H", N+1)]]














layer = [[("H", N+1)],[("CZ", (N+1, 4))],
        gatelayer("Rx", N; (θ=π*g)),
        [("XX+YY_couple", coupling_sequence[1][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[1])],
        [("XX+YY_couple", coupling_sequence[2][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[2])],
        [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
        [("CZ", (N+1, 4))],
        [("H", N+1)]]

layer = [[("H", N+1)],[("CZ", (N+1, 4))],
        gatelayer("Rx", N; (θ=π*g)),
        [("YY_couple", coupling_sequence[1][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[1])],
        [("YY_couple", coupling_sequence[2][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[2])],
        [("XX_couple", coupling_sequence[1][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[1])],
        [("XX_couple", coupling_sequence[2][i], (ϕ=-π/2,)) for i=1:length(coupling_sequence[2])],
        [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
        [("CZ", (N+1, 4))],
        [("H", N+1)]]
