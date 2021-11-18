using PastaQ
using ITensors
using Plots
using Statistics
using Distributions
using FFTW


# Arrange ZZ(:rand) couplings using CNOT-Rz-CNOT sequence
function ZZ_layer(N::Int)
    qarray = PastaQ.lineararray(N)
    first_sublattice_CNOTs = gatelayer("CX", qarray[1])
    second_sublattice_CNOTs = gatelayer("CX", qarray[2])
    randϕ = rand(Uniform(-π*1.5,-π*0.5), N) # set disordered angles for Ising coupling
    first_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 1:2:N]
    second_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 2:2:N]
    circuit = [[first_sublattice_CNOTs], [first_sublattice_Rzs], [first_sublattice_CNOTs], [second_sublattice_CNOTs], [second_sublattice_Rzs], [second_sublattice_CNOTs]]
    return circuit
end

Z_field(N::Int) = [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]

function sg_evolve(N, nsteps)
    t_vec = Vector{Float64}();
    χ_sg_final = Vector{Float64}(); # final average spin glass for each g
    g_vec = collect(0.7:0.02:1.0)

    println("Evolution running")

    for g=0.7:0.02:1.0 # varying g

        println("Running with g = $g")

        X_flip(N::Int) = gatelayer("Rx", 4; (θ=π*g)) # varied g for X_flip

        av_χ_sg = Vector{Float64}(); # average spin glass values for each run

        for r=1:10 # number of runs to average on
            println("   Run $r")

            χ_sg_vec = Vector{Float64}();

            #ψ = initial_bitstring
            ψ = productstate(N)

            for i=1:nsteps

                append!(t_vec, i)

                # applies TC unitary
                ψ = runcircuit(ψ, X_flip(N))
                normalize!(ψ)
                ψ = runcircuit(ψ, ZZ_layer(N))
                normalize!(ψ)
                ψ = runcircuit(ψ, Z_field(N))
                normalize!(ψ)

                zzcorr = correlation_matrix(ψ, "Z", "Z"; site_range=2:N-1) # correlation function excludes edge qubits

                χ_sg = 0

                # calculates χˢᵍ order parameter as in google paper
                for i=1:size(zzcorr)[1], j=1:size(zzcorr)[2]
                    if i == j
                        χ_sg += 0
                    else
                        χ_sg += 1/(N-2)*(zzcorr[i,j])^2
                    end
                end

                append!(χ_sg_vec, abs(χ_sg))
            end

            χ_sg_mean = mean(χ_sg_vec[50:60]) # average between t=50 and 5=60
            println(χ_sg_mean)
            append!(av_χ_sg, χ_sg_mean)
        end
        run_mean = mean(av_χ_sg)
        append!(χ_sg_final, run_mean) # takes average of all 10 runs for each g
    end
    return g_vec, χ_sg_final
end

N = 20 # no. of qubits
nsteps = 60 # no. of evolutions

initial_bitstring = randomstate(N) # generate random state

g_vec, χ_sg_final = sg_evolve(N, nsteps);
plot(g_vec, χ_sg_final, xaxis="g", yaxis="Spin Glass Order Parameter, χˢᵍ", label="N=$N", legend=true)
