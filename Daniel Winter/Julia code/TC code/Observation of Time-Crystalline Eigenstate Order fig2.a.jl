using Yao
using Random
using Plots
using Distributions
theme(:dao)

N=11
g=0.97
Zstring = chain(N, prod([put(N, i=>Z) for i=1:N]))
ZZpairs = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])
Xrotstring = chain(N, prod([put(N, i=>Rx(g)) for i=1:N]))
H =  Xrotstring + ZZpairs + Zstring
U(Jt::Float64) = time_evolve(H, Jt, tol=1e-5, check_hermicity=false)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

macro Name(arg)
   string(arg)
end

protected = true

function Mz_evolve(nsteps::Int64)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = rand_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            ψ |> U(0.01)
        else
            ψ |> U(0.01)
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
    mkpath(string(dir,"\\Graphs\\","\\",name,"\\","\\",q, " qubits\\","\\Protected ",protected,"\\","\\",step," steps"))
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",q, " qubits\\","\\Protected ",protected,"\\","\\",step," steps\\",strlabel,".png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

stp =100
i=0.000001

println("Successfully started "*string(stp)*" steps\n")

plot_name = @Name Obs_TC_fig2a

t_vec, Mz_vec = Mz_evolve(stp, i)

title1 = "" *string(plot_name)* " plot for "* string(N) * " Qubits for " * string(stp)* " cycles"
plt=Plots.plot(t_vec, Mz_vec, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, stp), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
Plots.title!(title1)
label=title1;

save_plot(plot_name, N, stp, plt, label) # Saves the plots to github
display(plt)

println("Successfully finished "*string(stp)*" steps\n")
