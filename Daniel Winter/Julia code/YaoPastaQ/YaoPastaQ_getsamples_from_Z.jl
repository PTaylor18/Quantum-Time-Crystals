using ITensors
using PastaQ
using Random
using HDF5
using Statistics
using Distributions
using Plots
theme(:dao)

Random.seed!(1234)

N = 20
depth = 2
nshots = 100
circuit = randomcircuit(N, depth)
# 1. Generation of measurement data on the quantum states
# at the output of a circuit. Each data-point is a projetive
# measurement in an arbitrary local basis. The default local basis
# is `["X","Y","Z"]`.
# a) Unitary circuit
# Returns output state as MPS
println("Generate samples from random projective measurements of the state U|0,0,…>:")
#data, ψ = getsamples(circuit, nshots; local_basis=["X", "Y", "Z"])
data, ψ = getsamples(circuit, nshots; local_basis=["Z"])
resetqubits!(ψ)

function getshots(circ, nshots, j)
    data, ψ = getsamples(circ, nshots; local_basis=["Z"])
    measurement_list = [];
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    results = Vector{Float64}();
    circuit = Vector[]
    qarray = PastaQ.lineararray(N)

    for i=1:nshots
        results = data[i, :]
        #push!(circuit, [("X", i)])
        #push!(circuit, [("CX", qarray[1])])
        for n=1:j
            measured_bit = parse(Float64, replace("Z", results[n]))
            measured_spin = 2 * (measured_bit - 0.5)
            append!(measurement_list, measured_spin)
        end
        append!(Mz_vec, mean(measurement_list))
        append!(t_vec, i)
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

plot_name = @Name getsamples

for q = N # Number of qubits for each
    for step = nshots
        figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/DTC " * string(q)* " qubits/" * string(step) * " steps/"
        # 1 after variable names denote they're local variables in the for loop
        # And here is where the file path is defined for each iteration.
        title1 = plot_name* " plot for "* string(q) * " Qubits for " * string(step)* " cycles"
        t_vec1, Mz_vec1 = getshots(circuit, nshots, N)
        plt=Plots.plot(t_vec1, Mz_vec1, linetype=:steppre, xlabel = "Time / period of deriving field ", xlims = (0, nshots), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
        Plots.title!(title1)
        #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

        label=title1;
        display(plt)
        save_plot(plot_name, q, step, plt, label) # Saves the plots to github
        println("Successfully finished "*string(step)*" steps\n")
    end
    println("Successfully finished "*string(q)*" qubits\n")
end
