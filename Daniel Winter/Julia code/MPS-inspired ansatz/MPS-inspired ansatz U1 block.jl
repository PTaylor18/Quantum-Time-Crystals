using Yao
using Random
using Plots
using FFTW
theme(:dao)

X1(N::Int) = chain(N, put(N, 1=>X))

CNOTodd(N::Int) = chain(N, prod([chain(N, cnot(i,i+1)) for i=1:N-1]))
CNOTend(N::Int) = chain(N, prod(chain(N, cnot(N,1))))

Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / (N)

RZZSwap(N::Int) = sum([chain(N, put(N, i=>Rz(0.1))*put(N,i+1=>Rz(0.1)))*chain(N, swap(N,i,i+1)) for i=1:N-1]) #*put(swap(N,i,i+1))
RZZend(N::Int) = sum([chain(N, put(N, 1=>Rz(0.1))*put(N, N=>Rz(0.1)))])
Swapend(N::Int) = chain(N, swap(N,1,N))

function Mz_evolve(N::Int, nsteps::Int64, Jt)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if isodd(i) == true
            ψ |> X1(N) |> RZZSwap(N) |> RZZend(N) |> Swapend(N)
        else
            ψ |> RZZSwap(N) |> RZZend(N) |> Swapend(N)

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
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

plot_name = @Name MPS-inspired-U1

for q = 5 # Number of qubits for each
    for step = 3
        figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/" *plot_name* " " *string(q)* " qubits/" * string(step) * " steps/"
        # 1 after variable names denote they're local variables in the for loop
        # And here is where the file path is defined for each iteration.
        for i = 0
            title1 = "" *string(plot_name)* " plot for "* string(q) * " Qubits for Jt = " * string(i/10)
            t_vec1, Mz_vec1 = Mz_evolve(q, step, i/10)
            plt=Plots.plot(t_vec1, Mz_vec1, linetype=:steppre, xlabel = "Time/ period of driving field", xlims = (0, step), ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
            Plots.title!(title1)
            #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

            label=title1;

            save_plot(plot_name, q, step, plt, label) # Saves the plots to github
            display(plt)

            tit1 = "Fourier transfrom of "*string(plot_name)* " plot for "* string(q) * " Qubits for " * string(step)* " cycles for Jt = " * string(i/10)
            Mz_vec_fft1 = fft(Mz_vec1)
            fourier=Plots.plot(abs.(Mz_vec_fft1), linetype=:steppre, xlabel="Frequecy", xlims = (0, step), ylabel="Intensity", legend = false)
            Plots.title!(tit1)
            #display(fourier)
            lab="Fourier transfrom of "*string(plot_name)* " plot for "* string(q) * " Qubits for " * string(step)* " cycles for Jt = " * string(i/10);

            #save_plot(plot_name, q, step, fourier, lab) # Saves the plots to github
        end
        println("Successfully finished "*string(step)*" steps\n")
    end
    println("Successfully finished "*string(q)*" qubits\n")
end

println("Successfully finished")
