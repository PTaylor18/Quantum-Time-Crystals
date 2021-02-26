using Yao
using Random
using Plots

Xstr(N::Int) = chain(N, repeat(X, 1:N))
Ystr(N::Int) = chain(N, repeat(Y, 1:N))
Zstr(N::Int) = chain(N, repeat(Z, 1:N))
CNOTodd(N::Int) = chain(N, prod([chain(N, cnot(i,i+1)) for i=1:2:N]))
CNOTeven(N::Int) = chain(N, prod([chain(N, cnot(i,(i+1)%N)) for i=2:2:N]))
Mzodd(N::Int) = sum([put(N, i => Z) for i = 1:2:N]) / (N/2)
Mzeven(N::Int) = sum([put(N, i => Z) for i = 2:2:N]) / (N/2)
Xstrodd(N::Int) = chain(N, repeat(X, 1:2:N))

RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.1)) for i=1:N]))
#Hzz(N::Int) = sum([-1 * put(N, i+1 => Z) * put(N, i => Z) for i = 1:N-1])
#Uzz(N::Int, Jt::Float64) = time_evolve(H_ZZ(N), Jt, tol=1e-5, check_hermicity=true)

function Mz_evolve(N::Int, nsteps::Int64, Jt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N) |> Xstrodd(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            ψ |> CNOTodd(N) |> CNOTeven(N) |> RXstr(N)
        else
            ψ |> CNOTodd(N) |> CNOTeven(N)  |> RXstr(N)
        end
        append!(Mz_vec, expect(Mzodd(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(6, 36, 1.)
plot(t_vec, Mz_vec, linetype=:steppre)
