using Yao
using Random
using Plots
theme(:dao)

Xstr(N::Int) = chain(N, repeat(X, 1:N))
RZstr(N::Int) = chain(N, prod([put(N, i+1=>Rz(2*pi)) for i=1:N-1]))
RZstr5(N::Int) = chain(N, put(N, 5=>Rz(2*pi)))
RZstr6(N::Int) = chain(N, put(N, 6=>Rz(2*pi)))

CNOT(N::Int) = chain(N, prod([chain(N, cnot(i,i+1)) for i=1:N-1]))
CNOT16(N::Int) = chain(N, chain(N, cnot(1,6)))
CNOT25(N::Int) = chain(N, chain(N, cnot(2,5)))

Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / (N)

one_d = false

function Mz_evolve(N::Int, nsteps::Int64)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if one_d == true
            ψ |> Xstr(N) |> CNOT(N) |> RZstr(N) |> CNOT(N)
        else
            ψ |> Xstr(N) |> CNOT(N) |> RZstr(N) |> CNOT(N) |> CNOT16(N) |> RZstr6(N) |> CNOT16(N) |> CNOT25(N) |> RZstr5(N) |> CNOT25(N)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

macro Name(arg)
   string(arg)
end

function save_plot(name, q, step, plot, label)
    cd(dirname(@__FILE__));
    dir = pwd();
    println("Saving plots to the directory in $dir","\\Graphs")
    strlabel = string(label);
    mkpath(string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\"))
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\",strlabel,".png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

if one_d ==true
    plot_name = @Name one_d_frame
else
    plot_name = @Name two_d_frame
end

for q = 9 # Number of qubits for each
    for step = 200:200:2000
        figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/" *plot_name* " " *string(q)* " qubits/" * string(step) * " steps/"
        # 1 after variable names denote they're local variables in the for loop
        # And here is where the file path is defined for each iteration.
        title1 = "" *string(plot_name)* " plot for "* string(q) * " Qubits"
        t_vec1, Mz_vec1 = @time Mz_evolve(q, step)
        plt=Plots.plot(t_vec1, Mz_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
        Plots.title!(title1)
        #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

        label=title1;

        save_plot(plot_name, q, step, plt, label) # Saves the plots to github
        display(plt)
        println("Successfully finished "*string(step)*" steps\n")
    end
    println("Successfully finished "*string(q)*" qubits\n")
end
