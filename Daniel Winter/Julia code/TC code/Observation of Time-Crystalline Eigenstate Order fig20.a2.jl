using Yao
using Random
using Plots
using Distributions
theme(:dao)


g=0.97

RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(π*g)) for i=1:N]))
ZZpairs(N) = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])
RZstr(N::Int) = chain(N, prod([put(N, i=>Rz(rand(Uniform(-π, π)))) for i=1:N]))

Hzz(N::Int) = sum([-1 * rand(Uniform(-π*1.5, -π*0.5)) * put(N, i+1 => Z) * put(N, i => Z) for i = 1:N-1])
Uzz(N::Int, Jt::Float64) = time_evolve(Hzz(N), Jt, tol=1e-5, check_hermicity=true)

Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

# magnetization on centre qubit
#centre_qubit = N/2
Mz_centre(N::Int) = chain(N, put(N, 3 => Z))

macro Name(arg)
   string(arg)
end

protected = true

function Mz_evolve(N::Int, nsteps::Int64, Jt)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            ψ |> RXstr(N) |> Uzz(N, Jt) |> RZstr(N)
        else
            ψ |> RXstr(N) |> RZstr(N)
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

plot_name = @Name Obs_TC_fig2a

for q = 11 # Number of qubits for each
    for step = 1000
        println("Successfully started " * string(step) *" Steps \n")
        figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/DTC " * string(q)* " qubits/" * string(step) * " steps/"
        # 1 after variable names denote they're local variables in the for loop
        # And here is where the file path is defined for each iteration.
        for i = 0
            title1 = "Plot for "* string(q) * " Qubits for " * string(step)* " cycles for Jt = " * string(i/10)
            t_vec1, Mz_vec1 = Mz_evolve(q, step, i/10)
            plt=Plots.plot(t_vec1[1:step], Mz_vec1[1:step], linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
            Plots.title!(title1)
            #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

            label=title1;
            display(plt)
            save_plot(plot_name, q, step, plt, label) # Saves the plots to github
        end
        println("Successfully finished "*string(step)*" steps\n")
    end
    println("Successfully finished "*string(q)*" qubits\n")
end
