using ITensors
using PastaQ
using Random
using Statistics

println("Generate samples from random projective measurements of the state U|0,0,…>:")
function getmagn(circ, nshots, N)
    data, ψ = getsamples(circ, nshots; local_basis=["Z"])
    Mz = 0.;
    for j=1:N
        results = data[:, j]
        measurement_list = [];
        for n=1:nshots
            measured_bit = parse(Float64, replace("Z", results[n]))
            measured_spin = 2 * (measured_bit - 0.5)
            append!(measurement_list, measured_spin)
        end
        Mz += mean(measurement_list)
    end
    return Mz
end

Random.seed!(1234)

N = 4
depth = 20
nshots = 100
flip(N::Int) = gatelayer("X", N)

#circuit = randomcircuit(N, depth)
circuit1 = flip(N)
measured_shots = getmagn(circuit1, nshots, N)

circuit2 = [flip(N), flip(N)]
measured_shots = getmagn(circuit2, nshots, N)
