using ITensors
using PastaQ
using Random
using HDF5

Random.seed!(1234)

N = 4
depth = 2
nshots = 10
circuit = randomcircuit(N, depth)

# 1. Generation of measurement data on the quantum states
# at the output of a circuit. Each data-point is a projetive
# measurement in an arbitrary local basis. The default local basis
# is `["X","Y","Z"]`.
# a) Unitary circuit
# Returns output state as MPS
println("Generate samples from random projective measurements of the state U|0,0,…>:")
#data, ψ = getsamples(circuit, nshots; local_basis=["X", "Y", "Z"])
data, ψ = getsamples(circuit, nshots; local_basis=["Z"])

function getshots(circ, nshots, j)
    data, ψ = getsamples(circ, nshots; local_basis=["Z"])
    results = data[:, j]
    measurement_list = [];
    for n=1:nshots
        measured_bit = parse(Float64, replace("Z", results[n]))
        measured_spin = 2 * (measured_bit - 0.5)
        append!(measurement_list, measured_spin)
    end
    return measurement_list
end
measured_shots = getshots(circuit, nshots, 1)
println(mean(measured_shots))
