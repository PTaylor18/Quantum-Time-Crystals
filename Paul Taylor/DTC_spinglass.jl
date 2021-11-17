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

initial_bitstring = randomstate(4)

function Mz_evolve_sg(N, nsteps)
    t_vec = Vector{Float64}();
    av_χ_sg = Vector{Float64}();
    g_vec = collect(0.7:0.02:1.0)

    for g=0.7:0.02:1.0 # varying g

        println("Running with g = $g")
        X_flip(N::Int) = gatelayer("Rx", 4; (θ=π*g))

        χ_sg_vec = Vector{Float64}();

        #ψ = initial_bitstring
        ψ = productstate(N)
        for i=1:nsteps

            append!(t_vec, i)

            ψ = runcircuit(ψ, X_flip(N))
            normalize!(ψ)
            ψ = runcircuit(ψ, ZZ_layer(N))
            normalize!(ψ)
            ψ = runcircuit(ψ, Z_field(N))
            normalize!(ψ)

            zzcorr = correlation_matrix(ψ, "Z", "Z"; site_range=2:N-1) # correlation function excludes edge qubits

            χ_sg = 0

            for i=1:size(zzcorr)[1], j=1:size(zzcorr)[2]
                if i == j
                    χ_sg += 0
                else
                    χ_sg += 1/(N-2)*(zzcorr[i,j])^2
                end
            end
            append!(χ_sg_vec, abs(χ_sg))
        end
        χ_sg_mean = mean(χ_sg_vec[50:60])
        append!(av_χ_sg, χ_sg_mean)
    end
    return g_vec, av_χ_sg
end

nsteps = 60

g_vec, av_χ_sg = Mz_evolve_sg(4, nsteps);
plot(g_vec, av_χ_sg, xaxis="g", yaxis="Spin Glass Order Parameter, χˢᵍ", legend=false)


plot(t_vec[1:nsteps], χ_sg_vec[1:nsteps], xaxis="Time T", yaxis="Spin Glass Order Parameter, χˢᵍ", legend=false)


χ_sg_mean = mean(χ_sg_vec[50:60])
append!(av_χ_sg, χ_sg_mean)

g_vec = collect(0.7:0.02:1.0)
g_vec = collect(0.7:0.1:1.0)



χᴳ = 0

for i=1:size(zzcorr)[1], j=1:size(zzcorr)[2]
    if i == j
        χᴳ += 0
    else
        χᴳ += zzcorr[i,j]
    end
end

χᴳ

size(zzcorr)[1]



zzcorr[1,2]
zzcorr[1,3]
zzcorr[1,4]
zzcorr[2,4]
zzcorr[1,1]
