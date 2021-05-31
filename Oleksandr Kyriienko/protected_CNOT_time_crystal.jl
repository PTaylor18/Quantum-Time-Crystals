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

RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.05)) for i=1:N]))
RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(rand(Uniform(0.05,0.15)))) for i=1:N]))


using Distributions
RZstr(N::Int) = chain(N, prod([put(N, i=>Rz(rand(Uniform(0.5,2.5)))) for i=1:N]))


#Hz(N::Int) = rand(Uniform(0.5,2.5)) * sum([put(N, i => Z) for i = 1:N]) + 0.05 * sum([put(N, i => X) for i = 1:N])
#RZstr(N::Int) = time_evolve(Hz(N), 1., tol=1e-5, check_hermicity=true)


function Mz_evolve(N::Int, nsteps::Int64, Jt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    #ψ = zero_state(N) |> Xstrodd(N)
    ψ = zero_state(N) |> chain(N, put(1 => X))
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            ψ |> CNOTodd(N) |> CNOTeven(N) |> RZstr(N) #|> RXstr(N)
        else
            ψ |> CNOTodd(N) |> CNOTeven(N)  |> RXstr(N)
        end
        append!(Mz_vec, expect(Mzodd(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(10, 500, 1.)
plot(t_vec, Mz_vec, linetype=:steppre)
