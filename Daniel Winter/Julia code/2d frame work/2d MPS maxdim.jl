using PastaQ
using ITensors
using Random
using Plots
using Statistics
using Distributions
using CurveFit # used for curve fitting

theme(:dao)

one_d = true # Is the system 1d
qmin = 4
qmax = 20

nq = Vector{Float64}() # number of qubits
t1 = Vector{Float64}() # time for 1d
t2 = Vector{Float64}() # time for 2d

# Define custom function to measure an observable, in this
# case a Pauli operator on `site`
function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

import PastaQ: gate

# Ising (ZZ) coupling gate
function gate(::GateName"ZZ"; ϕ::Number)
  return [
    exp(-im*ϕ/2) 0 0 0
    0 exp(im*ϕ/2) 0 0
    0 0 exp(im*ϕ/2) 0
    0 0 0 exp(-im*ϕ/2)
  ]
end

gate(::GateName"ZZ_couple"; kwargs...) = gate("ZZ"; kwargs...)

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
coupling_sequence13 = [(1,3)]
coupling_sequence24 = [(2,4)]
coupling_sequence2d = [(1,3),(2,4)]

function coupling_seq(N)
    sequence  = []
    for i=1:2:N-1
      push!(sequence, (i,i+1))
    end
    return sequence
end

function coupling_seq13(N)
    sequence13  = []
    for i=1
      push!(sequence13, (1,3))
    end
    return sequence13
end

function coupling_seq24(N)
    sequence24  = []
    for i=1
      push!(sequence24, (2,4))
    end
    return sequence24
end

function coupling_seq2d(N)
    sequence2d  = []
    for i=1:2:N÷2
      push!(sequence2d, (i,N÷2 + i))
    end
    return sequence2d
end

σz(ψ::MPS) = [measure_pauli(ψ, j, "Z") for j in 1:length(ψ)]

function Mz_evolve(N, nsteps)
    t_vec = Vector{Float64}();
    link_vec = Vector{Float64}();
    maxlink_vec = Vector{Float64}();
    coupling_sequence = coupling_seq(N);

    circuit = Vector[];

    layer = [gatelayer("X", N),
        [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)]
        ]

    for i=0:nsteps
        append!(t_vec, i)
        push!(circuit, layer)
    end

    # define the Circuit observer
    obs = Observer([
      "χs" => linkdims,      # bond dimension at each bond
      "χmax" => maxlinkdim,  # maximum bond dimension
      "σᶻ" => σz,
    ])
    ψ = runcircuit(circuit; (observer!)=obs)
    link_res = results(obs, "χs")
    maxlink_res = results(obs, "χmax")
    for i=1:length(link_res)
        append!(link_vec, mean(link_res[i]))
    end
    for i=1:length(maxlink_res)
        append!(maxlink_vec, mean(maxlink_res[i]))
    end
    return t_vec, link_vec, maxlink_vec
end

function Mz_evolve2d(N, nsteps)
    t_vec = Vector{Float64}();
    link_vec = Vector{Float64}();
    maxlink_vec = Vector{Float64}();
    coupling_sequence = coupling_seq(N);

    circuit = Vector[];

    #coupling_sequence13 = coupling_seq13(N);
    #coupling_sequence24 = coupling_seq24(N);
    coupling_sequence2d = coupling_seq2d(N);

    layer = [gatelayer("X", N),
    [("ZZ_couple", coupling_sequence[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence)],
    [("ZZ_couple", coupling_sequence2d[i], (ϕ=rand(Uniform(-π*1.5,-π*0.5)),)) for i=1:length(coupling_sequence2d)],
    ]

    for i=0:nsteps
        append!(t_vec, i)
        push!(circuit, layer)
    end

    # define the Circuit observer
    obs = Observer([
      "χs" => linkdims,      # bond dimension at each bond
      "χmax" => maxlinkdim,  # maximum bond dimension
      "σᶻ" => σz,
    ])
    ψ = runcircuit(circuit; (observer!)=obs)
    link_res = results(obs, "χs")
    maxlink_res = results(obs, "χmax")
    for i=1:length(link_res)
        append!(link_vec, mean(link_res[i]))
    end
    for i=1:length(maxlink_res)
        append!(maxlink_vec, mean(maxlink_res[i]))
    end
    return t_vec, link_vec, maxlink_vec
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
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\",strlabel," ",step," steps.png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

function save_comp_plot(plot, label)
    cd(dirname(@__FILE__));
    dir = pwd();
    println("Saving plots to the directory in $dir","\\Graphs")
    strlabel = string(label);
    mkpath(string(dir,"\\Graphs\\"))
    pngfilename = string(dir,"\\Graphs\\","Time comparison.png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

for runs = 1
    for q = qmin:1:qmax # Number of qubits for each
        for step = 50
            one_d = true

            if one_d == true
                plot_name = @Name one_d_frame
            else
                plot_name = @Name two_d_frame
            end
            figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/" *plot_name* " " *string(q)* " qubits/" * string(step) * " steps/"
            # 1 after variable names denote they're local variables in the for loop
            # And here is where the file path is defined for each iteration.
            title1 = "" *string(plot_name)* " maximum bond dimension plot for "* string(q) * " Qubits"
            t_vec1, link_vec1, maxlink_vec1 = Mz_evolve(q, step)
            plt=Plots.plot(t_vec1, maxlink_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Maximum bond dimension", legend = false)
            Plots.title!(title1)
            #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

            label=title1;
            display(plt)
            save_plot(plot_name, q, step, plt, label) # Saves the plots to github
            println("Successfully finished "*string(step)*" steps\n")
        end

        for step = 50
            one_d = false

            if one_d == true
                plot_name = @Name one_d_frame
            else
                plot_name = @Name two_d_frame
            end
            figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/" *plot_name* " " *string(q)* " qubits/" * string(step) * " steps/"
            # 1 after variable names denote they're local variables in the for loop
            # And here is where the file path is defined for each iteration.
            title1 = "" *string(plot_name)* " maximum bond dimension plot for "* string(q) * " Qubits"
            t_vec1, maxlink_vec1 = Mz_evolve2d(q, step)
            plt=Plots.plot(t_vec1, maxlink_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Maximum bond dimension", legend = false)
            Plots.title!(title1)
            #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

            label=title1;
            save_plot(plot_name, q, step, plt, label) # Saves the plots to github
            display(plt)
            println("Successfully finished "*string(step)*" steps\n")
        end

        println("Successfully finished "*string(q)*" qubits\n")

    end
end
