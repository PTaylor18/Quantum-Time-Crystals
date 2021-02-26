using Yao
using Random
using Plots

function Hising(N::Int, J::Float64)
    Hising = chain(N); # empty Hamiltonian to start
    for i=1:N-1
        Hising += chain(N, J * put(N, i=>Z) * put(N, (i+1)=>Z))
    end
    return Hising
end
HIsing = Hising(4, 1.)
mat(HIsing)

circuit = chain(4, cnot(1,2))
push!(circuit, chain(4, cnot(2,3)))
println(circuit)

function cnot_layer(N)
    circ = chain(N)
    for i=1:N-1
        push!(circ, chain(N, cnot(i, (i+1))))
    end
    circ
end
U=cnot_layer(4)

state = rand_state(4)
state |> U
statevec(state)

expect(HIsing, state) |> real


N = 4
mz = sum([put(N, i => Z) for i = 1:2:N])
expect(mz, state)
