using Yao
using Plots
using Compose
using YaoPlots

Xstring(N) = chain(N, prod([put(N, i=>X) for i=1:N]))
ZZpairs(N) = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])

RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.5)) for i=1:N]))

U(N::Int, Jt::Float64) = time_evolve(0.5 * Xstring(N) + 0.5 * ZZpairs(N), Jt, tol=1e-5, check_hermicity=true)
U_minus(N::Int, Jt::Float64) = time_evolve(-1*(0.5 * Xstring(N) + 0.5 * ZZpairs(N)), Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int64) = sum([put(N, i => Z) for i = 1:N]) / N

Corellation_Mz(N::Int, Jt::Float64) = chain(N, chain(N, put(N, i => Z) for i = 1:N),
    chain(N,time_evolve(0.5 * Xstring(N) + 0.5 * ZZpairs(N), Jt, tol=1e-5, check_hermicity=true)),
    chain(N, put(N, i => Z) for i = 1:N),
    chain(time_evolve(-1*(0.5 * Xstring(N) + 0.5 * ZZpairs(N)), Jt, tol=1e-5, check_hermicity=true)))

vizcircuit(Corellation_Mz(6,0.00001))

function Ghz_state(N::Int)
    ψ = chain(N, put(N, 1 => H))
    for i=1:N-1
        push!(ψ, chain(N, cnot(i, (i+1))))
    end
    ψ
end

vizcircuit(Ghz_state(3))

zero_state(3) |> Ghz_state(3)

statevec(zero_state(3) |> Ghz_state(3))

function Mz_evolve(N, nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    #ψ = zero_state(N)
    ψ = rand_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U(N, i * deltaJt) |> RXstr(N)
        else
            ψ |> U(N, i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

function Mz_evolve_Ghz(N::Int, nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N) |> Ghz_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U(N, i * deltaJt) |> RXstr(N)
        else
            ψ |> U(N, i * deltaJt)
        end
        append!(Mz_vec, real(expect(Corellation_Mz(N, π/4), ψ)))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(6,10000, 0.00001)
t_vec, Mz_vec = Mz_evolve_Ghz(6, 10000, 0.00001)
plot(t_vec[9000:10000], Mz_vec[9000:10000])
xlabel!("Time T")
ylabel!("Correlation Function")

#correlation function not working correctly atm

function zz_evolve(N::Int, angle::Float64)
    cirq = chain(N)
    for i=1:N-1
        # implement ZZ(θ) unitaries
        push!(cirq, chain(N, cnot(i,i+1)))
        push!(cirq, chain(N, put((i+1)=>Rz(2*angle))))
        push!(cirq, chain(N, cnot(i,i+1)))
    end
    return cirq
end

function xstr_evolve(N::Int, angle::Float64)
    cirq = chain(N)
    # implement XX..X(θ) string
    push!(cirq, repeat(H, 1:N))
    push!(cirq, chain(N, cnot(i,i+1) for i=1:N-1))
    push!(cirq, chain(N, put(N=>Rz(2*angle))))
    push!(cirq, chain(N, cnot(i,i+1) for i=1:N-1)')
    push!(cirq, repeat(H, 1:N))
    return cirq
end

function fast_Mz_evolve_Ghz(N, nsteps, deltaJt)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N) |> Ghz_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        ψ |> xstr_evolve(N, 0.5*i * deltaJt) |> zz_evolve(N, 0.5*i * deltaJt)
        append!(Mz_vec, real(expect(Corellation_Mz(N, i* deltaJt), ψ)))
    end
    return t_vec, Mz_vec
end

t_vec_fast, Mz_vec_fast = fast_Mz_evolve_Ghz(6, 10000, 0.00001)
plot(t_vec_fast[9000:10000], Mz_vec_fast[9000:10000])
xlabel!("Time T")
ylabel!("Correlation Function")

#paper on pauli matrices
#cuyao
