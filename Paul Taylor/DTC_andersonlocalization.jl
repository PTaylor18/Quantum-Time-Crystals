using PastaQ
using ITensors
using Plots
using Statistics
using Distributions
using FFTW

X_flip(N::Int) = gatelayer("Rx", 4; (θ=π*0.97))

Y_flip(N::Int) = gatelayer("Ry", 4; (θ=π*1.))

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

function ZZ_layer2(N::Int)
    qarray = PastaQ.lineararray(N)
    first_sublattice_CNOTs = gatelayer("CX", qarray[1])
    second_sublattice_CNOTs = gatelayer("CX", qarray[2])
    randϕ = rand(Uniform(-π*1.5,-π*0.5), N) # set disordered angles for Ising coupling
    first_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 1:2:N]
    second_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 2:2:N]
    circuit = [[first_sublattice_CNOTs], [first_sublattice_Rzs], [first_sublattice_CNOTs]]
    return circuit
end

ZZ_layer2(4)

function XXYY_layer(N::Int)
    qarray =  PastaQ.lineararray(N)
    first_lattice_iSwaps = gatelayer("iSwap", qarray[1])
    second_lattice_iSwaps= gatelayer("iSwap", qarray[2])
    circuit = [[first_lattice_iSwaps], [second_lattice_iSwaps]]
end

function XXYY_layer2(N::Int)
    qarray =  PastaQ.lineararray(N)
    first_lattice_iSwaps = gatelayer("iSwap", qarray[1])
    second_lattice_iSwaps= gatelayer("iSwap", qarray[2])
    circuit = [[second_lattice_iSwaps]]
end

XXYY_layer2(4)

Z_field(N::Int) = [("Rz", i, (ϕ=rand(Uniform(-π,π)),)) for i=1:N]

Y_field(N::Int) = [("Ry", i, (θ=rand(Uniform(-π,π)),)) for i=1:N]

XXYY_layer(4)

PastaQ.lineararray(4)

# function gate(::GateName"XX+YY"; ϕ::Number)
#   return [
#     (cos(-ϕ/2)^2 + sin(-ϕ/2)^2) 0 0 0
#     0 (-sin(-ϕ/2)^2 + cos(-ϕ/4)*cos(-ϕ/2)) im*(cos(-ϕ/2)*sin(-ϕ/2) - cos(-ϕ/4)*sin(-ϕ/2)) 0
#     0 (-2*im*cos(-ϕ/2)*sin(-ϕ/2)) (-sin(-ϕ/2)^2 + cos(-ϕ/2)^2) 0
#     0 0 0 (cos(-ϕ/2)^2 + sin(-ϕ/2)^2)
#   ]
# end

function gate(::GateName"XX+YY"; ϕ::Number)
  return [
    (cos(ϕ/2)^2 + sin(ϕ/2)^2) 0 0 0
    0 (cos(ϕ/2)^2 - sin(ϕ/2)*cos(ϕ/2)) (-2*im*sin(ϕ/2)*cos(ϕ/2)) 0
    0 (-2*im*cos(ϕ/2)*sin(ϕ/2)) (-sin(ϕ/2)^2 + cos(ϕ/2)^2) 0
    0 0 0 (cos(ϕ/2)^2 + sin(ϕ/2)^2)
  ]
end

gate(::GateName"XX+YY_couple"; kwargs...) = gate("XX+YY"; kwargs...)

coupling_sequence = coupling_seq(4)

test = [("XX+YY_couple", (1,2), (ϕ=rand(Uniform(-π*1.5,-π*0.5)),))]

test = [[("XX+YY_couple", coupling_sequence_test[1][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence_test[1])],
        [("XX+YY_couple", coupling_sequence_test[2][i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence_test[2])]]

function XXYY_layer_cirq(N::Int)

    coupling_sequence = PastaQ.lineararray(N)

    # circuit = [
    # [("XX+YY_couple", coupling_sequence[1][i], (ϕ=rand(Uniform(π*0.5,π)),)) for i=1:length(coupling_sequence[1])],
    # [("XX+YY_couple", coupling_sequence[2][i], (ϕ=rand(Uniform(π*0.5,π)),)) for i=1:length(coupling_sequence[2])]
    # ]

    circuit = [
    [("XX+YY_couple", coupling_sequence[1][i], (ϕ=-π,)) for i=1:length(coupling_sequence[1])],
    [("XX+YY_couple", coupling_sequence[2][i], (ϕ=-π,)) for i=1:length(coupling_sequence[2])]
    ]
end

XXYY_layer_cirq(4)

function neel_state(N)

    circuit = Vector[]

    for i=2:2:N
        one = [("X", i)]
        push!(circuit, one)
    end
    return circuit
end

neel_state(4)


function XXYY_evolve(N, nsteps)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();

    nshots = 50
    ψ = productstate(N)

    for i=1:nsteps
        println("Step: $i")
        append!(t_vec, i)

        ψ = runcircuit(ψ, X_flip(N))
        normalize!(ψ)
        ψ = runcircuit(ψ, XXYY_layer_cirq(N))
        normalize!(ψ)
        ψ = runcircuit(ψ, Z_field(N))
        normalize!(ψ)


        # ψ = runcircuit(ψ, X_flip(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, ZZ_layer(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, XXYY_layer(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, XXYY_layer_cirq(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, Z_field(N))
        # normalize!(ψ)

        # ψ = runcircuit(ψ, X_flip(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, ZZ_layer(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, XXYY_layer2(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, ZZ_layer2(N))
        # normalize!(ψ)
        # ψ = runcircuit(ψ, Z_field(N))
        # normalize!(ψ)

        data = getsamples(ψ, nshots; local_basis=["Z"])
        #data, ρ = getsamples(circuit, nshots; local_basis=["Z"])

        res = Vector{Float64}();

        for x in data
            append!(res, 2 * (x.second - 0.5))
        end

        append!(Mz_vec, mean(res))
    end
    return t_vec, Mz_vec
end

N = 4 # No. of qubits
nsteps = 100

t_vec, Mz_vec = XXYY_evolve(N, nsteps);
plot(t_vec[1:nsteps], Mz_vec[1:nsteps], linetype=:steppre, xaxis="Time T", yaxis="Magnetisation Mz", legend=false)

# plot fourier transform
Mz_vec_fft = fft(Mz_vec)
freq_vec = (t_vec.-1)./nsteps
plot(freq_vec, (abs.(Mz_vec_fft)), xaxis="Frequency (1/T)", yaxis="Fourier Spectrum", legend=false)
