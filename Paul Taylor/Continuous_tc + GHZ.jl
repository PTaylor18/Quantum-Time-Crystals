using Yao
using Compose
using YaoPlots

N = 6
Xstring = chain(N, prod([put(N, i=>X) for i=1:N]))
ZZpairs = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])

U(Jt::Float64) = time_evolve(0.5 * Xstring + 0.5 * ZZpairs, Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N
Mz_0(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

H_gates(N::Int) = chain(N, repeat(H, 1:N))
CNOTodd(N::Int) = chain(N, prod([chain(N, cnot(i,i+1)) for i=1:2:N]))
CNOTeven(N::Int) = chain(N, prod([chain(N, cnot(i,(i+1)%N)) for i=2:2:N-1]))

Ghz_state(N) = chain(N, repeat(H, 1:N),
    prod([chain(N, cnot(i,i+1)) for i=1:2:N]),
    prod([chain(N, cnot(i,(i+1)%N)) for i=2:2:N-1]),
    prod([chain(N, cnot(i,i+1)) for i=1:2:N]),
    prod([chain(N, cnot(i,(i+1)%N)) for i=2:2:N-1]),
    repeat(H, 1:N))

vizcircuit(Ghz_state(6))

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

function Mz_evolve_Ghz(nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            #ψ |> Ghz_state(N) |> U(i * deltaJt)
            ψ |> H_gates(N) |> CNOTodd(N) |> CNOTeven(N) |> U(i * deltaJt) |>
            CNOTodd(N) |> CNOTeven(N) |> H_gates(N)

        else
            ψ |> Ghz_state(N) |> U(i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

#Correlation(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N


using Plots
t_vec, Mz_vec = Mz_evolve(10000, 0.00001)
t_vec, Mz_vec = Mz_evolve_Ghz(10000, 0.00001)
plot(t_vec[9000:10000], Mz_vec[9000:10000])

expect
