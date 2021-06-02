using Yao
using Plots
using YaoPlots
using FFTW


Xstring(N) = chain(N, prod([put(N, i=>X) for i=1:N]))
Ystring(N) = chain(N, prod([put(N, i=>Y) for i=1:N]))
Zstring(N) = chain(N, prod([put(N, i=>Z) for i=1:N]))

ZZpairs(N) = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])
XXpairs(N) = sum([chain(N, put(N, i=>X)*put(N, i+1=>X)) for i=1:N-1])
YYpairs(N) = sum([chain(N, put(N, i=>Y)*put(N, i+1=>Y)) for i=1:N-1])

XX_YYpairs(N) = sum([chain(N, put(N, i=>X)*put(N, i+1=>X)) + put(N, i=>Y)*put(N, i+1=>Y) for i=1:N-1])

U_XXYY(N::Int, Jt::Float64) = time_evolve(0.5 * Xstring(N) + 0.5 * XX_YYpairs(N) + 0.5 * ZZpairs(N), Jt, tol=1e-5, check_hermicity=true)

RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.1)) for i=1:N]))

U(N::Int, Jt::Float64) = time_evolve(0.5 * Xstring(N) + 0.5 * ZZpairs(N), Jt, tol=1e-5, check_hermicity=true)
U_minus(N::Int, Jt::Float64) = time_evolve(-1*(0.5 * Xstring(N) + 0.5 * ZZpairs(N)), Jt, tol=1e-5, check_hermicity=true)

U_XYZstr(N::Int, Jt::Float64) = time_evolve(0.5 * Xstring(N) + 0.5 * Ystring(N)
    + 0.5 * Zstring(N) + 0.5 * ZZpairs(N), Jt, tol=1e-5, check_hermicity=true)

U_2(N::Int, Jt::Float64) = time_evolve(0.5 * Xstring(N) + 0.5 * ZZpairs(N) + 0.5 * XXpairs(N) + 0.5 * YYpairs(N), Jt, tol=1e-5, check_hermicity=true)

Mz(N::Int64) = sum([put(N, i => Z) for i = 1:N]) / N

vizcircuit(chain(6, put(6, i => Z) for i = 1:6))

Correlation_Mz(N::Int, Jt::Float64) = chain(N, chain(N, sum([put(N, i => Z) for i = 1:N]) / N),
    chain(N, time_evolve(0.5 * Xstring(N) + 0.5 * ZZpairs(N), Jt, tol=1e-5, check_hermicity=true)),
    chain(N, sum([put(N, i => Z) for i = 1:N]) / N),
    chain(N, time_evolve(-1*(0.5 * Xstring(N) + 0.5 * ZZpairs(N)), Jt, tol=1e-5, check_hermicity=true)))

vizcircuit(Correlation_Mz(6,0.00001))

function Ghz_state(N::Int)
    #Generates a Ghz state for N qubits
    ψ = chain(N, put(N, 1 => H))
    for i=1:N-1
        push!(ψ, chain(N, cnot(i, (i+1))))
    end
    ψ
end

vizcircuit(Ghz_state(6))
zero_state(6) |> Ghz_state(6)
statevec(zero_state(6) |> Ghz_state(6))

using Random

RX_1qubit(N) = chain(N, put(N, 6 => Rx(0.1))) #Small x rotation on 1 qubit that acts as a perturbation

H_state(N::Int) = chain(N, repeat(H, 1:N))

function Mz_evolve(N, nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    #ψ = rand_state(N)
    #ψ = zero_state(N) |> H_state(N)
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

t_vec, Mz_vec = Mz_evolve(6,10000, 0.00001)
plot(t_vec[1:5000], Mz_vec[1:5000])
xlabel!("Time T")
ylabel!("Average Magnetization, Mz")

function Mz_evolve_Ghz(N::Int, nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N) |> Ghz_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U(N, i * deltaJt) |> RX_1qubit(N)
        else
            ψ |> U(N, i * deltaJt)
        end
        append!(Mz_vec, real(expect(Correlation_Mz(N, i * 0.00001), ψ)))
    end
    return t_vec, Mz_vec
end


t_vec, Mz_vec = Mz_evolve_Ghz(6, 10000, 0.00001)
plot(t_vec[1:10000], Mz_vec[1:10000])
xlabel!("Time T")
ylabel!("Correlation Function <Mz(t)Mz(0)>")

#correlation function not working correctly atm

Mz_vec_fft = fft(Mz_vec)
plot(abs.(Mz_vec_fft[1:100]))


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
Mz(N::Int64) = sum([put(N, i => Z) for i = 1:N]) / N

test(N) = chain(N, [put(N,i => X) for i = 1:N-3])
vizcircuit(test(6))

Xstring_half(N) = chain(N, prod([put(N, i=>X) for i=1:3]))
ZZpairs_half(N) = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:5])

U_half(N::Int, Jt::Float64) = time_evolve(0.5 * Xstring_half(N) + 0.5 * ZZpairs_half(N), Jt, tol=1e-5, check_hermicity=true)

function Mz_evolve_half(N, nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    #ψ = rand_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U_half(N, i * deltaJt) |> RXstr(N)
        else
            ψ |> U(N, i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve_half(6,10000, 0.00001)
plot(t_vec[1:5000], Mz_vec[1:5000])
xlabel!("Time T")
ylabel!("Average Magnetization, Mz")
