using Yao
# using Compose
# using YaoPlots
using Plots
using Distributions

N = 6 # number of qubits

RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(π*0.97)) for i=1:N]))
ZZpairs(N) = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])
RZstr(N::Int) = chain(N, prod([put(N, i=>Rz(rand(Uniform(-π, π)))) for i=1:N]))

Hzz(N::Int) = sum([-1 * rand(Uniform(-π*1.5, -π*0.5)) * put(N, i+1 => Z) * put(N, i => Z) for i = 1:N-1])
Uzz(N::Int, Jt::Float64) = time_evolve(Hzz(N), Jt, tol=1e-5, check_hermicity=true)

Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

# magnetization on centre qubit
centre_qubit = N/2
Mz_centre(N::Int) = chain(N, put(N, 3 => Z))
function Mz_evolve(N::Int, nsteps::Int64, Jt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    #ψ = product_state(bit"01001")
    for i = 0:nsteps
        append!(t_vec,i)
        if protected
            ψ |> RXstr(N) |> Uzz(N, Jt) |> RZstr(N)
            # ψ |> U_google(N, Jt)
        else
            ψ |> XstrR(N) |> RXstr(N)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(N, 10000, 1.)
plot(t_vec[1:1000], Mz_vec[1:1000], linetype=:steppre,  xaxis=("Time T"), yaxis=("Magnetisation Mz"), legend=false)
