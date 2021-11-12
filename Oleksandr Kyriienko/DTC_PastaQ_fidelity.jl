using ITensors, PastaQ
using PastaQ: MPS, MPO
using Distributions
using Plots
using PastaQ: resetqubits!

function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

# define the total magnetization measurement for pure states
Mz_pure(N::Int, ψ::MPS) = sum([measure_pauli(ψ, j, "Z") for j in 1:N]) # /N

# To define magentization as an MPO we use ITensors.jl `AutoMPO()` function
function Mz_mixed(N::Int, rho::MPO)
    sites = siteinds("Qubit", N)
    ampo = AutoMPO()
    for j in 1:N
      ampo .+= 1., "Z", j
    end
    Mz_MPO = MPO(ampo, sites)
    # Measure the total magnetization
    Mz = real(tr(rho * Mz_MPO))
    return Mz
end

# Apply the X gates on each qubit
flip(N::Int) = gatelayer("X", N) # == [("X", j) for j in 1:N]

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

# evolve initial state in time using the stroboscopic DTC sequence
function Mz_evolve(N::Int, nsteps::Int, γX::Float64, γZZ::Float64)
    # Initialize the MPS state ψ = |0,0,..,0⟩
    #ψ = productstate(N);
    # Write magnetization from each step (starting with 0)
    Mz_list = [1.];
    #ρ = dot(productstate(N), productstate(N)'); # this variable stores MPO if there is some noise
    #ρ = copy(ψ);
    #ρ = productstate(N);
    ψ0 = randomstate(N) #productstate(N)
    #ψ1 = runcircuit(ψ0, flip(N))
    if (γX > 1E-8) && (γZZ > 1E-8)
        # Initialize the MPS state ρ = |0,0,..,0><0,0,..0|
        sites = siteinds("Qubit", N);
        #ρ = LPDO(productstate(N));
        #ρ0 = LPDO(productstate(N));
        ρ = MPO(randomstate(sites; χ = 10, ξ = 15));
        ρ0 = MPO(randomstate(sites; χ = 10, ξ = 15));
        #resetqubits!(ψ0)
        resetqubits!(ρ)
        resetqubits!(ρ0)
        ρ1 = runcircuit(ρ0, flip(N))
        #fidelity
        #ρ = productoperator(N);
        for n in 1:nsteps
            ρ = runcircuit(ρ, flip(N); noise = ("amplitude_damping", (γ = γX,))) # , apply_dag=true
            #results = measure(ρ, "Z"); println(results)
            normalize!(ρ) # normalize the density operator
            ρ = runcircuit(ρ, ZZ_layer(N); noise = ("amplitude_damping", (γ = γZZ,))) # , apply_dag=true
            normalize!(ρ) # normalize the density operator
            #println(PastaQ.array(ρ))
            #measurement = Mz_mixed(N, ρ)
            isodd(n) ? (measurement = -1. * fidelity(ρ, ρ1)) : (measurement = 1. * fidelity(ρ, ρ0))
            println(measurement)
            append!(Mz_list, measurement)
        end
    else
        ψ = randomstate(N) #productstate(N); # Initialize the MPS state ψ = |0,0,..,0⟩
        for n in 1:nsteps
            ψ = runcircuit(ψ, flip(N))
            ψ = runcircuit(ψ, ZZ_layer(N))
            append!(Mz_list, Mz_pure(N, ψ))
        end
    end

    return Mz_list
end

N = 10 # number of qubits
nsteps = 100 # number of Trotter steps
γX =0.001 # single-qubit error rate
γZZ = 0.0001 # two-qubit error rate
Mz_DTC = Mz_evolve(N, nsteps, γX, γZZ)
plt = plot(collect(0:nsteps), Mz_DTC, linetype=:steppre, xaxis=("discrete time"), yaxis=("Mz(t)"), label="N=$N", legend=true)

println("Finished")
