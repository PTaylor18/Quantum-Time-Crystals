using PastaQ
using ITensors
using Random
using Plots
theme(:dao)

# Set the random seed so the results are the same
# each time it is run
Random.seed!(1234)

ψ = qubits(4)

# Apply the X gate on qubit N
gX(N::Int) = ("X", N)
# Apply the Y gate on qubit N
gY(N::Int) = ("Y", N)
# Apply the Z gate on qubit N
gZ(N::Int) = ("Z", N)
# Apply the Controlled-X on qubits N
gCXodd(N::Int) = ("CX", (N, N + 1))
gCXeven(N::Int) = ("CX", (N, (N+1)%N))
# Apply the Rotation of θ around X on qubits N
gRx(N::Int) = ("Rx", N, (θ=0.1,))

gCNOTodd(N::Int) = ("CNOT", (N, N + 1))
gCNOTeven(N::Int) = ("CNOT",(N, (N+1)%N))


function Mz_evolve(N::Int, nsteps::Int64)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();

    # Initialize the MPS state ψ = |0,0,0,...,0⟩
    ψ = qubits(N)
    for a=1:2:N
        ψ = runcircuit(ψ,gX(a))
    end
    for i = 0:nsteps
        append!(t_vec, i)
        for j=1:2:N
            ψ = runcircuit(ψ,gCNOTodd(j))
        end
        for k=2:2:N
            ψ = runcircuit(ψ,gCNOTeven(k))
        end
        for l=1:N
            ψ = runcircuit(ψ,gRx(l))
        end
        append!(Mz_vec, sum(getsamples(ψ,1))/N)
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
    mkpath(string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\","\\",step," steps"))
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

plot_name = @Name QuantumGate

for q = 10 # Number of qubits for each
    for step = 1000
        figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/DTC " * string(q)* " qubits/" * string(step) * " steps/"
        # 1 after variable names denote they're local variables in the for loop
        # And here is where the file path is defined for each iteration.
        title1 = plot_name* " plot for "* string(q) * " Qubits for " * string(step)* " cycles"
        t_vec1, Mz_vec1 = Mz_evolve(q, step)
        plt=Plots.plot(t_vec1, Mz_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
        Plots.title!(title1)
        #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

        label=title1;
        display(plt)
        save_plot(plot_name, q, step, plt, label) # Saves the plots to github
        println("Successfully finished "*string(step)*" steps\n")
    end
    println("Successfully finished "*string(q)*" qubits\n")
end
