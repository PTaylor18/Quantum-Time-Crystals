#Quantum Packages
using Yao
using Random
using Plots

#Menu Packages
using TerminalMenus

options = ["1. DTC", "2. CNOT",]

menu = RadioMenu(options, pagesize=3)

choice = request("What type of crystal are you testing? :\n", menu)

if choice != -1
    println("Using the ", options[choice], "\n")


    Xstr(N::Int) = chain(N, repeat(X, 1:N))
    Ystr(N::Int) = chain(N, repeat(Y, 1:N))
    Zstr(N::Int) = chain(N, repeat(Z, 1:N))
    CNOTodd(N::Int) = chain(N, prod([chain(N, cnot(i,i+1)) for i=1:2:N]))
    CNOTeven(N::Int) = chain(N, prod([chain(N, cnot(i,(i+1)%N)) for i=2:2:N]))
    Mzodd(N::Int) = sum([put(N, i => Z) for i = 1:2:N]) / (N/2)
    Mzeven(N::Int) = sum([put(N, i => Z) for i = 2:2:N]) / (N/2)
    Xstrodd(N::Int) = chain(N, repeat(X, 1:2:N))

    RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.1)) for i=1:N]))
    #Hzz(N::Int) = sum([-1 * put(N, i+1 => Z) * put(N, i => Z) for i = 1:N-1])
    #Uzz(N::Int, Jt::Float64) = time_evolve(H_ZZ(N), Jt, tol=1e-5, check_hermicity=true)
    if choice == "1"

        function Mz_evolve(N::Int, nsteps::Int64, Jt)
            protected = true
            t_vec = Vector{Float64}();
            Mz_vec = Vector{Float64}();
            ψ = zero_state(N) |> Xstrodd(N)
            for i = 0:nsteps
                append!(t_vec, i)
                if protected
                    ψ |> CNOTodd(N) |> CNOTeven(N) |> RXstr(N)
                else
                    ψ |> CNOTodd(N) |> CNOTeven(N)  |> RXstr(N)
                end
                append!(Mz_vec, expect(Mzodd(N), ψ))
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
            mkpath(string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps"))
            pngfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".png")
            #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
            #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
            savefig(plot, pngfilename)
            #savefig(plot, pdffilename)
            #savefig(plot, epsfilename)
        end

        plot_name = @Name CNOT_TC

        for q = 2:2:6 # Number of qubits for each
            for step = 10:5:15
                figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/" *plot_name* " " *string(q)* " qubits/" * string(step) * " steps/"
                # 1 after variable names denote they're local variables in the for loop
                # And here is where the file path is defined for each iteration.
                for i = 0:20
                    title1 = "" *string(plot_name)* " plot for "* string(q) * " Qubits for " * string(step)* " cycles for Jt = " * string(i/10)
                    t_vec1, Mz_vec1 = Mz_evolve(q, step, i/10)
                    plt=Plots.plot(t_vec1, Mz_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
                    Plots.title!(title1)
                    #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

                    label=title1;

                    save_plot(plot_name, q, step, plt, label) # Saves the plots to github
                end
                println("Successfully finished "*string(step)*" steps\n")
            end
            println("Successfully finished "*string(q)*" qubits\n")
        end

        println("Successfully finished")

    else
        println("Fail!")
    end

else
    println("Menu canceled.")
end
