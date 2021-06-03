using Yao

#Listing 12: CUDA register
using CuYao
using CUDA

#constructs the |1010> state
#r = ArrayReg(bit"1010");

#transfer data to CUDA
#r = cu(ArrayReg(bit"1010"));

#Listing 13: instruct! and measure
#r = zero_state(4);
#instruct!(r, Val(:X), (1, ))
#samples = measure(r; nshots=3)
#[samples[1]...]


N = 6
Xstring = chain(N, prod([put(N, i=>X) for i=1:N]))
ZZpairs = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])

U(Jt::Float64) = time_evolve(0.5 * Xstring + 0.5 * ZZpairs, Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function Mz_evolve(nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U(i * deltaJt)
        else
            ψ |> U(i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

using Plots
t_vec, Mz_vec = Mz_evolve(10000, 0.00001)
plot(t_vec[9000:10000], Mz_vec[9000:10000])


# check if the operators commute
iscommute(Xstring, ZZpairs)


# implement the sequence of commuting ZZ evolutions as CNOT - Rz - CNOT
function zz_evolve(n::Int, angle::Float64)
    cirq = chain(n)
    for i=1:n-1
        # implement ZZ(θ) unitaries
        push!(cirq, chain(n, cnot(i,i+1)))
        push!(cirq, chain(n, put((i+1)=>Rz(2*angle))))
        push!(cirq, chain(n, cnot(i,i+1)))
    end
    return cirq
end

function xstr_evolve(n::Int, angle::Float64)
    cirq = chain(n)
    # implement XX..X(θ) string
    push!(cirq, repeat(H, 1:n))
    push!(cirq, chain(n, cnot(i,i+1) for i=1:n-1))
    push!(cirq, chain(n, put(n=>Rz(2*angle))))
    push!(cirq, chain(n, cnot(i,i+1) for i=1:n-1)')
    push!(cirq, repeat(H, 1:n))
    return cirq
end


# do some tests
using Random
Jt = rand()
psi = rand_state(N)
psi1 = copy(psi) |> zz_evolve(N, Jt)
psi2 = copy(psi) |> time_evolve(ZZpairs, Jt, tol=1e-10, check_hermicity=true)
fidelity(psi1, psi2) #fidelity is the measure of the distance between two quantum states.

Zstring = chain(N, prod([put(N, i=>Z) for i=1:N]))

Jt = rand()
psi = rand_state(N)
psi1 = copy(psi) |> xstr_evolve(N, Jt)
psi2 = copy(psi) |> time_evolve(Xstring, Jt, tol=1e-10, check_hermicity=true)
fidelity(psi1, psi2)


function fast_Mz_evolve(nsteps, deltaJt)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        ψ |> xstr_evolve(N, 0.5*i * deltaJt) |> zz_evolve(N, 0.5*i * deltaJt)
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec_fast, Mz_vec_fast = fast_Mz_evolve(10000, 0.00001)
plot(t_vec_fast[9000:10000], Mz_vec_fast[9000:10000])


function pauli_string_exp(n, dict, theta)
    cirq = chain(n)
    dict_keys = [key for key in keys(dict)]

    for key in dict_keys
        if dict[key] == "X"
            push!(cirq, put(n, key=>H)  )
        elseif dict[key] == "Y"
            push!(cirq, put(n, key=>S')  )
            push!(cirq, put(n, key=>H)  )
        end
    end

    for i_key in 2:length(dict_keys)
        push!(cirq, cnot(n, dict_keys[i_key - 1], dict_keys[i_key]))
    end

    push!(cirq, chain(n, put(dict_keys[end]=>Rz(theta))) )

    for i_key in length(dict_keys):-1:2
        push!(cirq, cnot(n, dict_keys[i_key - 1], dict_keys[i_key]))
    end

    for key in dict_keys
        if dict[key] == "X"
            push!(cirq, put(n, key=>H)  )
        elseif dict[key] == "Y"
            push!(cirq, put(n, key=>H)  )
            push!(cirq, put(n, key=>S)  )
        end
    end

    return cirq
end
d = Dict(1=>"X", 2=>"X", 3=>"X", 4=>"X", 5=>"X", 6=>"X");
d[2]
d_keys = [key for key in keys(d)]

pauli_string_exp(6, d, 0.1)

Jt = rand()
psi = rand_state(N)
psi1 = copy(psi) |> xstr_evolve(N, Jt)
#psi2 = copy(psi) |> time_evolve(Xstring, Jt, tol=1e-10, check_hermicity=true)
psi2 = copy(psi) |> pauli_string_exp(N, d, 2*Jt)
fidelity(psi1, psi2)
