using PastaQ

# Example 2: generation of measurement data

# Set parameters
N = 4                           # Number of qubits
depth = 4                       # Depth of random circuit
nshots = 1000                   # Number of measurements
gates = randomcircuit(N, depth) # Build gates


# 2a) Output state of a noiseless circuit. By default, each projective measurement
# is taken in basis randomly drawn from the the Pauli group. Also returns the output MPS.
data, ψ = getsamples(N, gates, nshots)

#  Note: the above is equivalent to:
# > bases = randombases(N, nshots; local_basis = ["X","Y","Z"])
# > ψ = runcircuit(N, gates)
# > data = getsamples(ψ, bases)

# 2b) Output state of a noisy circuit. Also returns the output MPO
data, ρ = getsamples(N, gates, nshots; noise = ("amplitude_damping", (γ = 0.01,)))

# 2c) Generate data for quantum process tomography, consisting of input states
# to a quantum channel, and the corresponding projective measurements
# at the output. By defaul, the states prepared at the inputs are selected from
# product states of eigenstates of Pauli operators, while measurements bases are
# sampled from the Pauli group.

# Unitary channel, returns the MPO unitary circuit
data, U = getsamples(N, gates, nshots; process=true)

# Noisy channel, returns the Choi matrix
data, Λ = getsamples(N, gates, nshots; process = true, noise = ("amplitude_damping", (γ = 0.01,)))
