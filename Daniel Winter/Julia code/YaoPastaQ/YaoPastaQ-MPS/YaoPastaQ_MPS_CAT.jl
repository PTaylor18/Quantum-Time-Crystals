using PastaQ
using ITensors
using Printf
using Distributions
using ITensors
import PastaQ: gate
using Plots
theme(:dao)


# Apply the H gate on qubit N
gH(N::Int) = ("H", N)
# Apply the CNOT gate on qubit N
gCNOTodd(N::Int) = ("CNOT", (N, N + 1))
gCNOTeven(N::Int) = ("CNOT",(N, (N+1)%N))


function gate(::GateName"H";)
  return [
    1/√2 1/√2
    1/√2 -1/√2
  ]
end

function gate(::GateName"CNOT";)
    return [
    1 0 0 0
    0 1 0 0
    0 0 0 1
    0 0 1 0
    ]
end

# Define custom function to measure an observable, in this
# case a Pauli operator on `site`
function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

# coupling indices for ZZ gates
coupling_sequence = [(1,2),(3,4)]

function coupling_seq(N)
    sequence  = []
    for i=1:2:N-1
      push!(sequence, (i,i+1))
    end
    return sequence
end

# create the TC Unitary
circuit = [gatelayer("H", 1),
        [("CNOT", coupling_sequence[i]) for i=1:length(coupling_sequence)]]

σx2(ψ::MPS) = measure_pauli(ψ, 2, "Z")
σz05(ψ::MPS) = measure_pauli(ψ, 1, "Z")
σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]

function Mz_evolve(N, nsteps)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    coupling_sequence = coupling_seq(N);
    circuit = Vector[];
    # first layer
    for j=1:N
        append!(t_vec, j)
        push!(circuit, gatelayer("H", 1))
    end
    layer = [[("CNOT", (i, i+1)) for i=1:N-1],
            [("CNOT", (k, k+1)) for k=N-1:-1:1]]
    for i=N:nsteps
        append!(t_vec, i)
        push!(circuit, layer)
    end

    # define the Circuit observer
    obs = Observer([
      "χs" => linkdims,      # bond dimension at each bond
      "χmax" => maxlinkdim,  # maximum bond dimension
      "σᶻ(0.5)" => σz05,        # pauli X on site .5
      "σᶻ" => σz,
    ])
    ψ = runcircuit(circuit; (observer!)=obs)
    Mz_res = results(obs, "σᶻ")
    for i=1:length(Mz_res)
        append!(Mz_vec, mean(Mz_res[i]))
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

plot_name = @Name NoisyCAT

for q = 2:10 # Number of qubits for each
    for step = 100
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
