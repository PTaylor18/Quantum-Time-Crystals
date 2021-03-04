using Yao
#using YaoPlots
using Compose
using Plots

Xstr(N::Int) = chain(N,repeat(X, 1:N))

Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function crz_layer(N)
    circ = chain(N)
    for i=1:N-1
        push!(circ, control(N, i,i+1=>Rz(0.5)))
    end
    circ
end
U=crz_layer(4)
#vizcircuit(U)

function crx_layer(N, perturbation)
    circ = chain(N)
    for i=1:N-1
        push!(circ, control(N, i,i+1=>Rx(perturbation)))
    end
    circ
end

function Mz_evolve(N::Int, nsteps::Int64, perturbation::Float64)
    protected = false
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec,i)
        if protected
            ψ |> Xstr(N) |> crz_layer(N)
        else
            ψ |> Xstr(N) |> crz_layer(N) |> crx_layer(N, perturbation)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(3, 60, 0.5)
plot(t_vec, Mz_vec, linetype=:steppre)
xlabel!("Time")
ylabel!("Magnetisation Mz")
