using PastaQ

# Example 1b: noisy quantum circuit

N = 4     # Number of qubits
depth = 4 # Depth of the quantum circuit
gates = randomcircuit(N, depth) # random circuit

# Run the circuit using an amplitude damping channel with decay rate `γ=0.01`.
# Returns the MPO for the mixed density operator `ρ = ε(|0,0,…⟩⟨0,0,̇…|), where
# `ε` is the quantum channel.
ρ = runcircuit(N, gates; noise = ("amplitude_damping", (γ = 0.01,)))

# Compute the Choi matrix of the channel
Λ = runcircuit(N, gates; process = true, noise = ("amplitude_damping", (γ = 0.01,)))
