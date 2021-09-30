using Yao
using Random
using Plots
using Distributions
theme(:dao)

N=6
Zstring = chain(N, prod([put(N, i=>Z) for i=1:N]))
ZZpairs = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])
Xrotstring = chain(N, prod([put(N, i=>Rx(0.05)) for i=1:N]))

U(Jt::Float64) = time_evolve(0.5 * Zstring + 0.5 * ZZpairs + 0.5 * Xrotstring, Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

macro Name(arg)
   string(arg)
end

plot_name = @Name TimeEvo_CNOT_rand_state

protected = true

function Mz_evolve(nsteps::Int64, deltaJt)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = rand_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            ψ |> U(i * deltaJt)
        else
            ψ |> U(i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec, ψ
end


function save_plot(name, q, step, plot, label)
    cd(dirname(@__FILE__));
    dir = pwd();
    println("Saving plots to the directory in $dir","\\Graphs")
    strlabel = string(label);
    mkpath(string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\","\\Protected ",protected,"\\","\\",step," steps"))
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\","\\Protected ",protected,"\\","\\",step," steps\\",strlabel,".png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

t_vec, Mz_vec = Mz_evolve(10000, 0.00001)
plot(t_vec[9000:10000], Mz_vec[9000:10000])
