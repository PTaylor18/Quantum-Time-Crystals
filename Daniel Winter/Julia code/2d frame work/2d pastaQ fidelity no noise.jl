using ITensors, PastaQ
using PastaQ: MPS, MPO
using Plots
using Distributions
using PastaQ: resetqubits!

using CurveFit # used for curve fitting
#using LsqFit # Causes issues with PastaQ

theme(:dao)

one_d = true # Is the system 1d
γX = 1E-2 # single-qubit error rate
γZZ = 1E-4 # two-qubit error rate
res = 1E-8 # The cut off

nq = Vector{Float64}() # number of qubits
t1 = Vector{Float64}() # time for 1d
t2 = Vector{Float64}() # time for 2d

function measure_pauli(ψ::MPS, site::Int, pauli::String)
  ψ = orthogonalize!(copy(ψ), site)
  ϕ = ψ[site]
  obs_op = gate(pauli, firstsiteind(ψ, site))
  T = noprime(ϕ * obs_op)
  return real((dag(T) * ϕ)[])
end

# define the total magnetization measurement for pure states
Mz_pure(N::Int, ψ::MPS) = sum([measure_pauli(ψ, j, "Z") for j in 1:N]) # /N

# To define magentization as an MPO we use ITensors.jl `AutoMPO()` function
function Mz_mixed(N::Int, rho::MPO)
    sites = siteinds("Qubit", N)
    ampo = AutoMPO()
    for j in 1:N
      ampo .+= 1., "Z", j
    end
    Mz_MPO = MPO(ampo, sites)
    # Measure the total magnetization
    Mz = real(tr(rho * Mz_MPO))
    return Mz
end

# Apply the X gates on each qubit
flip(N::Int) = gatelayer("X", N) # == [("X", j) for j in 1:N]

# Arrange ZZ(:rand) couplings using CNOT-Rz-CNOT sequence
function ZZ_layer(N::Int)
    qarray = PastaQ.lineararray(N)
    first_sublattice_CNOTs = gatelayer("CX", qarray[1])
    second_sublattice_CNOTs = gatelayer("CX", qarray[2])
    randϕ = rand(Uniform(-π*1.5,-π*0.5), N) # set disordered angles for Ising coupling
    first_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 1:2:N]
    second_sublattice_Rzs = [("Rz", j, (ϕ=randϕ[j],)) for j in 2:2:N]
    circuit = [[first_sublattice_CNOTs], [first_sublattice_Rzs], [first_sublattice_CNOTs], [second_sublattice_CNOTs], [second_sublattice_Rzs], [second_sublattice_CNOTs]]
    return circuit
end

# Arrange ZZ(:rand) couplings using CNOT-Rz-CNOT sequence for the 2 and 5 interaction
function ZZ_3layer(N::Int)
    qarray = PastaQ.lineararray(N)
    first_sublattice_CNOTs = gatelayer("CX", qarray[1])
    second_sublattice_CNOTs = gatelayer("CX", qarray[2])
    randϕ = rand(Uniform(-π*1.5,-π*0.5), N) # set disordered angles for Ising coupling
    first_sublattice_Rzs = [("Rz", 1, (ϕ=randϕ[1],))]
    second_sublattice_Rzs = [("Rz", 3, (ϕ=randϕ[3],))]
    circuit = [[first_sublattice_CNOTs], [first_sublattice_Rzs], [first_sublattice_CNOTs], [second_sublattice_CNOTs], [second_sublattice_Rzs], [second_sublattice_CNOTs]]
    return circuit
end

# Arrange ZZ(:rand) couplings using CNOT-Rz-CNOT sequence for the 1 and 6 interaction
function ZZ_4layer(N::Int)
    qarray = PastaQ.lineararray(N)
    first_sublattice_CNOTs = gatelayer("CX", qarray[1])
    second_sublattice_CNOTs = gatelayer("CX", qarray[2])
    randϕ = rand(Uniform(-π*1.5,-π*0.5), N) # set disordered angles for Ising coupling
    first_sublattice_Rzs = [("Rz", 2, (ϕ=randϕ[2],))]
    second_sublattice_Rzs = [("Rz", 4, (ϕ=randϕ[4],))]
    circuit = [[first_sublattice_CNOTs], [first_sublattice_Rzs], [first_sublattice_CNOTs], [second_sublattice_CNOTs], [second_sublattice_Rzs], [second_sublattice_CNOTs]]
    return circuit
end

# evolve initial state in time using the stroboscopic DTC sequence
function Mz_evolve(N::Int, nsteps::Int, γX::Float64, γZZ::Float64)
    # Initialize the MPS state ψ = |0,0,..,0⟩
    #ψ = productstate(N);
    # Write magnetization from each step (starting with 0)
    Mz_list = [1.];
    #ρ = dot(productstate(N), productstate(N)'); # this variable stores MPO if there is some noise
    #ρ = copy(ψ);
    #ρ = productstate(N);
    ψ0 = randomstate(N)
    #ψ1 = runcircuit(ψ0, flip(N))
    if (γX > 1E-8) && (γZZ > 1E-8)
        # Initialize the MPS state ρ = |0,0,..,0><0,0,..0|
        sites = siteinds("Qubit", N);
        #ρ = LPDO(productstate(N));
        #ρ0 = LPDO(productstate(N));
        ρ = MPO(randomstate(sites; χ = 10, ξ = 10, cutoff = res));
        ρ0 = MPO(randomstate(sites; χ = 10, ξ = 10, cutoff = res));
        #resetqubits!(ψ0)
        resetqubits!(ρ)
        resetqubits!(ρ0)
        ρ1 = runcircuit(ρ0, flip(N), cutoff = res)
        #fidelity
        #ρ = productoperator(N);
        for n in 1:nsteps
            if one_d == true
                ρ = runcircuit(ρ, flip(N), cutoff = res) # , apply_dag=true
                #results = measure(ρ, "Z"); println(results)
                normalize!(ρ) # normalize the density operator
                ρ = runcircuit(ρ, ZZ_layer(N), cutoff = res) # , apply_dag=true
                normalize!(ρ) # normalize the density operator
                #println(PastaQ.array(ρ))
                #measurement = Mz_mixed(N, ρ)
                isodd(n) ? (measurement = -1. * fidelity(ρ, ρ1)) : (measurement = 1. * fidelity(ρ, ρ0))
                append!(Mz_list, measurement)
            else
                ρ = runcircuit(ρ, flip(N), cutoff = res) # , apply_dag=true
                #results = measure(ρ, "Z"); println(results)
                normalize!(ρ) # normalize the density operator
                ρ = runcircuit(ρ, ZZ_layer(N), cutoff = res) # , apply_dag=true
                normalize!(ρ) # normalize the density operator
                ρ = runcircuit(ρ, ZZ_3layer(N), cutoff = res) # , apply_dag=true
                normalize!(ρ) # normalize the density operator
                ρ = runcircuit(ρ, ZZ_4layer(N), cutoff = res) # , apply_dag=true
                normalize!(ρ) # normalize the density operator
                #println(PastaQ.array(ρ))
                #measurement = Mz_mixed(N, ρ)
                isodd(n) ? (measurement = -1. * fidelity(ρ, ρ1)) : (measurement = 1. * fidelity(ρ, ρ0))
                append!(Mz_list, measurement)
            end
        end
    else
        ψ = randomstate(N); # Initialize the MPS state ψ = |0,0,..,0⟩
        for n in 1:nsteps
            if one_d == true
                ψ = runcircuit(ψ, flip(N), cutoff = res)
                ψ = runcircuit(ψ, ZZ_layer(N), cutoff = res)
                append!(Mz_list, Mz_pure(N, ψ))
            else
                ψ = runcircuit(ψ, flip(N), cutoff = res)
                ψ = runcircuit(ψ, ZZ_layer(N), cutoff = res)
                ψ = runcircuit(ψ, ZZ_3layer(N), cutoff = res)
                ψ = runcircuit(ψ, ZZ_4layer(N), cutoff = res)
                append!(Mz_list, Mz_pure(N, ψ))
            end
        end
    end

    return Mz_list
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

qmin = 4
qmax = 6

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
            title1 = "" *string(plot_name)* " plot for "* string(q) * " Qubits"
            Mz_DTC = Mz_evolve(q, step, γX, γZZ)
            append!(nq,q)
            append!(t1,@elapsed Mz_evolve(q, step, γX, γZZ))
            plt=Plots.plot(collect(0:step), Mz_DTC, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
            Plots.title!(title1)
            #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

            label=title1;

            save_plot(plot_name, q, step, plt, label) # Saves the plots to github
            display(plt)

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
            title1 = "" *string(plot_name)* " plot for "* string(q) * " Qubits"
            Mz_DTC = Mz_evolve(q, step, γX, γZZ)
            append!(t2,@elapsed Mz_evolve(q, step, γX, γZZ))
            plt=Plots.plot(collect(0:step), Mz_DTC, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
            Plots.title!(title1)
            #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

            label=title1;
            save_plot(plot_name, q, step, plt, label) # Saves the plots to github
            display(plt)
            println("Successfully finished "*string(step)*" steps\n")
        end

        println("Successfully finished "*string(q)*" qubits\n")

        title_c = "Time for computation"
        pltt=Plots.scatter(nq,t1,color="cyan4", xlabel = "Number of Qubits", ylabel = " \n"*"Time / seconds", legend = false)
        pltt=Plots.scatter!(nq,t2,color="orangered3")
        fit1 = curve_fit(ExpFit, nq, t1)
        fit2 = curve_fit(ExpFit, nq, t2)
        x=qmin:0.01:qmax
        pltt=Plots.plot!(x, fit1.(x),color="cyan4")
        pltt=Plots.plot!(x, fit2.(x),color="orangered3")
        Plots.title!(title_c)
        save_comp_plot(pltt, title_c)
        display(pltt)
    end
end

# Trying to make a better curve fit using LsqFut
# m(t, p) = p[1] + p[2] * exp.(p[3] * t + p[4])
# p0 = [0.5, 0.5, 0.5, 0.5]
# fitf = curve_fit(m, nq, t1, p0)

# model(t, p) = p[1] * exp.(-p[2] * t) +p[3]
# tdata = 0:0.1:20
# ydata = model(tdata, [1.0 2.0 3.0]) + 0.01*randn(length(tdata))
# p0 = [0.5, 0.5]

# fitff = curve_fit(model, tdata, ydata, p0)
# param = fitff.param

# pltff=Plots.plot(tdata, param[1]* exp.(-param[2] * tdata),color="orangered3")
