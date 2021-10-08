using PastaQ

# Example 1a: random quantum circuit

N = 4     # Number of qubits
depth = 4 # Depth of the circuit

# Generate a random quantum circuit built out of layers of single-qubit random
# rotations + CX gates, alternating between even and of odd layers.
gates = randomcircuit(N, depth)

@show gates

# Returns the MPS at the output of the quantum circuit: `|ψ⟩ = Û|0,0,…,0⟩`
ψ = runcircuit(N, gates)

# Generate the MPO for the unitary circuit:
U = runcircuit(N, gates; process=true)
