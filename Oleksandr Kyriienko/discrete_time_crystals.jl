using Yao
using Random
using Plots

Xstr(N::Int) = chain(N, prod([put(N, i=>X) for i=1:N]))
RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.3)) for i=1:N]))
#== or ==#
XstrR(N::Int) = chain(N, repeat(X, 1:N))
Hzz(N::Int) = sum([-1 * put(N, i+1 => Z) * put(N, i => Z) for i = 1:N-1])
Uzz(N::Int, Jt::Float64) = time_evolve(Hzz(N), Jt, tol=1e-5, check_hermicity=true) # It was H_ZZ but changed to Hzz as it wasn't defined but H_zz is
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function Mz_evolve(N::Int, nsteps::Int64, Jt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            ψ |> XstrR(N) |> Uzz(N, Jt) |> RXstr(N)
        else
            ψ |> XstrR(N) |> RXstr(N)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(4, 36, 1.)
plot(t_vec, Mz_vec, linetype=:steppre)


using YaoExtensions

c = variational_circuit(6, 6)

gatecount(c)
42+72
