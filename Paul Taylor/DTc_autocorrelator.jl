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

            randϕ = rand(Uniform(-0.75π, -0.25π), N + (-1)^(4%3)) # j disorder
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

#---

function Correlator_evolve(N::Int, g::Float64, j::Float64, σ::Float64,
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
            h_disorder_array = predefined_disorder_instance[2][k]
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
                randθ = h_disorder_array[i,:]

                # q controls which qubit is measured for auto-correlator measurement

                layer = [[("H", N+1)],[("CZ", (N+1, q))],
                        gatelayer("Rx", N; (θ=π*g)),
                        [("ZZ_couple", coupling_sequence[1][i], (ϕ=randϕ[i],)) for i=1:length(coupling_sequence[1])],
                        [("ZZ_couple", coupling_sequence[2][i], (ϕ=randϕ[length(coupling_sequence[1])+i],)) for i=1:length(coupling_sequence[2])],
                        [("Rz", i, (ϕ=randθ[i],)) for i=1:N],
                        [("CZ", (N+1, q))],
                        [("H", N+1)]]

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
r = 20
j = -0.15
σ = abs.(0.2*j)
σ = 0.2
nsteps = 60

initial_state = "random";

#---
# Simulations:

t_vec = collect(1:100)

predefined_disorder_instance = disorder_arrays(N, r, nsteps);
t_vec, correlator_array_pol, final_av_correlator_vec_pol, bit_string_pol = Correlator_evolve(N, g, j, σ, nsteps, r, true, "polarized");
t_vec, correlator_array_neel, final_av_correlator_vec_neel, bit_string_neel = Correlator_evolve(N, g, j, σ, nsteps, r, true, "neel");
t_vec, correlator_array_rand, final_av_correlator_vec_rand, bit_string_rand = Correlator_evolve(N, g, j, σ, nsteps, r, true, "random");

t_vec, correlator_array_08, final_av_correlator_vec_08, bit_string_08 = Correlator_evolve(N, 0.8, j, σ, nsteps, r, true, "polarized");
t_vec, correlator_array_94, final_av_correlator_vec_94, bit_string_94 = Correlator_evolve(N, 0.94, j, σ, nsteps, r, true, "polarized");
t_vec, correlator_array_g1, final_av_correlator_vec_g1, bit_string_g1 = Correlator_evolve(N, 1.0, j, σ, nsteps, r, true, "polarized");

begin
    plot(t_vec[1:nsteps], final_av_correlator_vec_08[1:nsteps],
        xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
        xtickfontsize=9,ytickfontsize=9, xguidefontsize=14, yguidefontsize=14,
        label=latexstring("g = 0.80"),
        color=:green, linewidth=4, linetype=:steppre, grid=false, legend=:topright, legendfontsize=11)
    plot!(t_vec[1:nsteps], final_av_correlator_vec_94[1:nsteps],
        linetype=:steppre, linewidth=2, label=latexstring("g = 0.94"), color=:red)
    plot!(t_vec[1:nsteps], final_av_correlator_vec_g1[1:nsteps],
        linetype=:steppre, label=latexstring("g = 1.00"), color=:blue)
end
savefig("autocorrelator.png")
# plot fourier transform
begin
    ac_pol_vec_fft = normalize!(fft(final_av_correlator_vec_08[10:nsteps]))
    ac_neel_vec_fft = normalize!(fft(final_av_correlator_vec_94[10:nsteps]))
    ac_rand_vec_fft = normalize!(fft(final_av_correlator_vec_g1[10:nsteps]))

    freq_vec = (t_vec[10:nsteps].-10)./(nsteps-9)
    plot(freq_vec, (abs.(ac_pol_vec_fft)), xaxis=L"\omega/\omega_D", yaxis=L"FFT",
    xtickfontsize=9,ytickfontsize=9, xguidefontsize=14, yguidefontsize=14,
    label=latexstring("g = 0.80"), color=:green, grid=false, legend=:topright, legendfontsize=11)
    plot!(freq_vec, (abs.(ac_neel_vec_fft)), label=latexstring("g = 0.94"), color=:red)
    plot!(freq_vec, (abs.(ac_rand_vec_fft)), label=latexstring("g = 1.00"), color=:blue )
    plot!(size=(530,400))
end
savefig("autocorrelatorfft.png")

t_vec[10:nsteps].-9
plot(t_vec[1:60], final_av_correlator_vec_pol[1:60],
    xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
    label="Polarized", xtickfontsize=9,ytickfontsize=9, xguidefontsize=14, yguidefontsize=14, dpi=200,
    marker=(:circle), markersize=5, color=:blueviolet, linewidth=1, grid=false, legend=:topright)

plot(t_vec[1:50], final_av_correlator_vec_pol[1:50],
    xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
    label=latexstring("N = $N, g = $g"), title=L"MBL \ DTC",
    color=:red3, linewidth=1, linetype=:steppre, grid=false, legend=:topright)

begin
    ac_pol_vec_fft = normalize!(fft(final_av_correlator_vec_pol[10:nsteps]))
    freq_vec = (t_vec[10:nsteps].-10)./(nsteps-10)
    plot(freq_vec, (abs.(ac_pol_vec_fft)), xaxis=L"\omega/\omega_D", yaxis=L"FFT",
        label=latexstring("g = $g"), color=:red3, grid=false, legend=:topright)
    plot!(size=(250,250))
end

plot(t_vec[1:30], final_av_correlator_vec_g1[1:30],
    xaxis=L"\textbf{t}", yaxis=L"\textbf{\overline{\langle\hat{Z}\:(0)\hat{Z}\:(t)\rangle}}",
    label=latexstring("N = $N, g = 1.0"), title=L"MBL \ DTC",
    color=:blueviolet, linewidth=1, linetype=:steppre, grid=false, legend=:topright)

begin
    ac_g1_vec_fft = normalize!(fft(final_av_correlator_vec_g1[10:nsteps]))
    freq_vec = (t_vec[10:nsteps].-10)./(100-10)
    plot(freq_vec, (abs.(ac_g1_vec_fft)), xaxis=L"\omega/\omega_D", yaxis=L"FFT",
        label=latexstring("g = 1.0"), color=:red3, grid=false, legend=:topright)
    plot!(size=(250,250))
end

plot!(t_vec[1:60], final_av_correlator_vec_neel[1:60],
    label="Néel", marker=(:circle), markersize=5, linewidth=1.5, color=:red3)

plot!(t_vec[1:60], final_av_correlator_vec_rand[1:60],
    label="Random", marker=(:circle), markersize=5, linewidth=2.5, color=:black)

savefig("initialstates.png")

# plot fourier transform
begin
    ac_pol_vec_fft = normalize!(fft(final_av_correlator_vec_pol))
    ac_neel_vec_fft = normalize!(fft(final_av_correlator_vec_neel))
    ac_rand_vec_fft = normalize!(fft(final_av_correlator_vec_rand))

    freq_vec = (t_vec.-1)./100
    plot(freq_vec, (abs.(ac_pol_vec_fft)), xaxis="Frequency (1/T)", yaxis="Fourier Spectrum",
    label="Polarized", color=:purple, grid=false, legend=:topright, dpi = 200)
    plot!(freq_vec, (abs.(ac_neel_vec_fft)), label="Néel", color=:red)
    plot!(freq_vec, (abs.(ac_rand_vec_fft)), label="Random", color=:black )
end
savefig("initialstatesfft.png")
d = Normal(j, abs.(0.2*j))
plot(x->pdf(d,x), xaxis = "x", yaxis="Probability Density", legend=false)

#--- Tests


correlator_array = zeros(1,nsteps)
av_auto_correlator_vec = Vector{Float64}();

coupling_sequence = PastaQ.lineararray(N)


Mz_auto_correlator_vec = Vector{Float64}();


N = 8
g = 0.97
nsteps = 100
begin
    random_bitstring = [parse(Int, x) for x in rand(("0", "1"), N)] # generates a random bistring
    bitstring_ancilla = append!(random_bitstring, 0) # add 0 for ancilla qubit
    random_bitstring_state = [bitstring_ancilla[n] == 1 ? ("X", n) : ("I", n) for n in 1:N+1] # generates bitstring circuit

    n = convert(Vector{Int64}, ones(N))
    neel_string = [isodd(i) == true ? 0 : 1 for i in 1:N]
    neel_ancilla = append!(neel_string, 0)
    neel_state = [neel_ancilla[j] == 1 ? ("X", j) : ("I", j) for j in 1:N+1]

    #ψ = runcircuit(neel_state)
    ψ = productstate(N+1)
    σz_auto_correlation(ψ::MPS) = measure_pauli(ψ, N+1, "Z") # x-axis projection of ancilla qubit prepared in |+X> (google paper)
    σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:(length(ψ)-1)]

    # layer = [gatelayer("Rx", N; (θ=(π*g))),
    #         [("ZZ_couple", coupling_sequence[1][i], (ϕ=0,)) for i=1:length(coupling_sequence[1])],
    #         [("ZZ_couple", coupling_sequence[2][i], (ϕ=0,)) for i=1:length(coupling_sequence[2])],
    #         [("Rz", i, (ϕ=rand(Uniform(0,π)),)) for i=1:N]]


    # layer = [[("H", N+1)],[("CZ", (N+1, 4))],
    #         gatelayer("Rx", N; (θ=π*g)),
    #         [("CZ", (N+1, 4))],
    #         [("H", N+1)]]

    obs = Observer([
      "σᶻ_auto-correlator" => σz_auto_correlation, # <ψ|Ẑ(0)Ẑ(t)|ψ> on selected qubit measured with ancilla qubit
      "σᶻ" => σz,
      ])

    circuit = Vector[];
    push!(circuit, neel_state)
    #push!(circuit, neel_state)

    Mz_auto_correlator_res = Vector{Float64}();

    for i in 1:100

        layer = [[("H", N+1)],[("CZ", (N+1, 4))],
                gatelayer("Rx", N; (θ=π*g)),
                [("ZZ_couple", coupling_sequence[1][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[1])],
                [("ZZ_couple", coupling_sequence[2][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[2])],
                [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
                [("CZ", (N+1, 4))],
                [("H", N+1)]]

        ψ = runcircuit(ψ, layer)

        append!(Mz_auto_correlator_res, PastaQ.measure(ψ, ("Z", N+1)))
        #push!(circuit, layer)
    end
    circuit

    #ψ = runcircuit(circuit; (observer!)=obs)
    println("Initial bitstring: ", prod(string.(random_bitstring)))
    #Mz_auto_correlator_res = results(obs, "σᶻ_auto-correlator")

    plot(collect(1:100), Mz_auto_correlator_res[1:100], linetype=:steppre)
end

ψ

PastaQ.measure(ψ, ("Z", N+1))

Mz_res = results(obs, "σᶻ")
mz_q_array = zeros(1,nsteps)

for q=1:N
    mz_q_vec = Vector{Float64}()
    for i=1:nsteps
        append!(mz_q_vec, Mz_res[i][q])
    end
    mz_q_array = vcat(mz_q_array, mz_q_vec')
end

mz_q_array = mz_q_array[Not(1),:]
begin
    q = 2 # select qubit
    plot(t_vec[1:nsteps], mz_q_array[q,:][1:nsteps],
        linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz",
        label="Qubit $q, g=$g", legend=:topright)
end



randomcircuit(N, 4; twoqubitgates="CZ", onequbitgates="Rn")
begin
    random_bitstring = [parse(Int, x) for x in rand(("0", "1"), N)] # generates a random bistring
    bitstring_ancilla = append!(random_bitstring, 0) # add 0 for ancilla qubit
    random_bitstring_state = [bitstring_ancilla[n] == 1 ? ("X", n) : ("I", n) for n in 1:N+1] # generates bitstring circuit

    n = convert(Vector{Int64}, ones(N))
    neel_string = [isodd(i) == true ? 0 : 1 for i in 1:N]
    neel_ancilla = append!(neel_string, 0)
    neel_state = [neel_ancilla[j] == 1 ? ("X", j) : ("I", j) for j in 1:N+1]

    σz_auto_correlation(ψ::MPS) = measure_pauli(ψ, N+1, "Z") # x-axis projection of ancilla qubit prepared in |+X> (google paper)
    σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:(length(ψ)-1)]

    circuit = Vector[];
    push!(circuit, neel_state)

    ψ = productstate(N+1)
end

for i in 1:1

    layer = [[("H", N+1)],[("CZ", (N+1, 4))],
            gatelayer("Rx", N; (θ=π*g)),
            [("ZZ_couple", coupling_sequence[1][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[1])],
            [("ZZ_couple", coupling_sequence[2][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence[2])],
            [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N],
            [("CZ", (N+1, 4))],
            [("H", N+1)]]

    #push!(circuit, layer)
end
circuit
ψ = runcircuit(circuit)
ψ = runcircuit(ψ, layer)
random_bitstring_state

PastaQ.measure(ψ, ("Z", N+1))

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

initial_state = "random"

ini_state, bit_string = initial_state_gen(N, initial_state)

bit_string

pol_vec = convert(Vector{Int64}, zeros(N))
bit_string = prod(string.(pol_vec))
