using PastaQ
using ITensors
using Printf
import PastaQ: gate


gate(::GateName"expτZ"; τ::Float64) = exp(τ * h * gate("Z"))
gate(::GateName"expτZZ"; τ::Float64) = exp(τ * Φ *kron(gate("Z"), gate("Z")))
gate(::GateName"expτXr"; τ::Float64, B::Float64) = exp(τ * pi * g * gate("X"))

N = 10    # Number of spins
h = -pi
Φ = -0.5*pi
g = 1
β = 1.0   # Inverse temperature
τ = 0.005 # Trotter step

# Depth of the circuit
depth = β ÷ τ

# Transverse field
z_layer = [("expτZ", j, (τ=τ, h=h)) for j in 1:N]
# Ising interactions
zz_layer = [("expτZZ", (j, j + 1), (τ=τ,)) for j in 1:(N - 1)]
# x rotated by πg
xr_layer = [("expτXr", j, (τ=τ, Φ=Φ)) for j in 1:N]

# Build the gate structure
circuit = []
for d in 1:depth
  append!(circuit, z_layer)
  append!(circuit, zz_layer)
  append!(circuit, xr_layer)
end

sites = siteinds("Qubit", N)
ampo = AutoMPO()
for j in 1:N
  # Transverse field
  ampo .+= -h, "Z", j
end
for j in 1:(N - 1)
  # Ising ZZ interactions
  ampo .+= -1, "Z", j, "Z", j + 1
end
for j in 1:N
  # Transverse field X
  ampo .+= -Φ, "X", j
end
# Generate Hamilotnian MPO
H = MPO(ampo, sites)

# Density-matrix renormalization group
dmrg_iter = 5      # DMRG steps
dmrg_cutoff = 1E-10   # Cutoff
Ψ0 = randomMPS(sites) # Initial state
sweeps = Sweeps(dmrg_iter)
maxdim!(sweeps, 10, 20, 30, 40, 50, 100)
cutoff!(sweeps, dmrg_cutoff)
# Run
println("Running DMRG to get ground state of transverse field Ising model:")
E, Ψ = dmrg(H, Ψ0, sweeps)
@printf("\nGround state energy:  %.8f  \n", E)
println("\n---------------------------------------\n")

#
# 2b. Run the imaginary-time circuit
#

β = 5.0 # Inverse temperature
Δ = 0.5  # Intermediate time-step
depth = Δ ÷ τ # Depth of the circuit
steps = β ÷ Δ # Total number of circuit application

# Initialize the density operator
ρ = PastaQ.identity_mpo(H)

println("Running imaginary time evolution to approximate the density matrix ρ = exp(-βH):")
for b in 1:steps
  # Run the circuit
  global ρ = runcircuit(ρ, circuit; cutoff=1E-12)

  # Normalize the density operatorr
  normalize!(ρ)

  # Measure the energy
  E_th = inner(ρ, H)

  @printf("β = %.1f : tr(ρH) = %.8f \n", (Δ * b), E_th)
end
