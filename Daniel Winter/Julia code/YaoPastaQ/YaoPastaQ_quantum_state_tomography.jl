using PastaQ

# Load the training data, as well as the target quantum state from file.
data, target = loadsamples("PATH_TO_DATAFILE.h5")
N = size(data, 2) # Number of qubits

# 1. Reconstruction with a variational wavefunction:
#
# Initialize a variational MPS with bond dimension χ = 10.
ψ0 = randomstate(N; χ = 10)

# Initialize stochastic gradient descent with learning rate η = 0.01
opt = SGD(η = 0.01)

# Run quantum state tomography
ψ = tomography(data, ψ0; optimizer = opt, target = target)

# 2. Reconstruction with a variational density matrix:
#
# Initialize a variational LPDO with bond dimension χ = 10 and Kraus dimension ξ = 2.
ρ0 = randomstate(N; mixed = true, χ = 10, ξ = 2)

# Run quantum state tomography
ρ = tomography(data, ρ0; optimizer = opt, target = target)
