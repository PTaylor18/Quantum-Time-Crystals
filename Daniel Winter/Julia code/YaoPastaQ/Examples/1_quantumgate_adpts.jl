using PastaQ
using ITensors
using Random
using Plots
theme(:dao)

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

# Define custom function to measure an observable, in this
# case a Pauli operator on `site`
function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

# pauli X on site 2
σx2(ψ::MPS) = measure_pauli(ψ, 2, "X")
σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]

# define the Circuit observer
obs = Observer([
  "χs" => linkdims,      # bond dimension at each bond
  "χmax" => maxlinkdim,  # maximum bond dimension
  "σˣ(2)" => σx2,        # pauli X on site 2
  "σᶻ" => σz,
])

function Mz_evolve(N::Int, nsteps::Int64)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();

    # Initialize the MPS state ψ = |0,0,0,...,0⟩
    circuit = qubits(N)
    # Initialize the MPS state ψ = |0,1,0,1,...,0,1>
    for a=1:2:N
        circuit = runcircuit(circuit,gX(a))
    end
    for i = 0:nsteps
        append!(t_vec, i)
        for j=1:2:N
            circuit = runcircuit(circuit,gCNOTodd(j))
        end
        for k=2:2:N
            circuit = runcircuit(circuit,gCNOTeven(k))
        end
        for l=1:N
            circuit = runcircuit(circuit,gRx(l))
        end
        ψ = runcircuit(circuit; (observer!)=obs)
        append!(Mz_vec, sum(sum(results(obs, "σᶻ"))))
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
