using Yao
using Random
using Plots


Xstr(N::Int) = chain(N, prod([put(N, i=>X) for i=1:N]))
RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.3)) for i=1:N]))
#== or ==#
XstrR(N::Int) = chain(N, repeat(X, 1:N))
Hzz(N::Int) = sum([-1 * put(N, i+1 => Z) * put(N, i => Z) for i = 1:N-1])
Uzz(N::Int, Jt::Float64) = time_evolve(Hzz(N), Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function Mz_evolve(N::Int, nsteps::Int64, Jt)
    protected = true # Can change this to run the alternative XstrR and Rxstr
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

for q = 2:5 # Number of qubits for each
    for step = 20:20:80
        figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/DTC " * string(q)* " qubits/" * string(step) * " steps/"
        for i = 0:50
            # 1 after variable names denote they're local variables in the for loop
            title1 = "DTC plot for "* string(q) * " Qubits for " * string(step)* " cycles for Jt = " * string(i/10)
            t_vec1, Mz_vec1 = Mz_evolve(q, step, i/10)
            Plots.plot(t_vec1, Mz_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
            Plots.title!(title1)
            Plots.savefig(figpath1*title1*".png")
        end
    end

end

# I have no idea what the relevance of this last bit is...
using YaoExtensions

c = variational_circuit(6, 6)

gatecount(c)
