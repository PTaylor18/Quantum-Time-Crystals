using Yao
using Random
using Plots
using FFTW
using Distributions
theme(:dao)

N=6
Xstring = chain(N, prod([put(N, i=>X) for i=1:N]))
ZZpairs = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])

U(Jt::Float64) = time_evolve(0.5 * Xstring + 0.5 * ZZpairs, Jt, tol=1e-16, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

macro Name(arg)
   string(arg)
end

plot_name = @Name Continuous_TC

protected = true

function Mz_evolve(N::Int, nsteps::Int64, deltaJt)
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U(i * deltaJt)
        else
            ψ |> U(i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end


function save_plot(name, N, step, plot, label)
    cd(dirname(@__FILE__));
    dir = pwd();
    println("Saving plots to the directory in $dir","\\Graphs")
    strlabel = string(label);
    mkpath(string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",N, " qubits\\","\\Protected ",protected,"\\","\\",step," steps"))
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",N, " qubits\\","\\Protected ",protected,"\\","\\",step," steps\\",strlabel,".png")
    #pdffilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".pdf")
    #epsfilename = string(dir,"\\Graphs\\","\\",name,"_",q, " qubits\\","\\",step," steps\\",strlabel,".eps")
    savefig(plot, pngfilename)
    #savefig(plot, pdffilename)
    #savefig(plot, epsfilename)
end

for j = N # Number of qubits for each
    for step = 10000
        figpath1 = "C:/Users/Daniel/OneDrive/Documents/Exeter Uni/Modules/Year 3/Project-Time crystals/Julia Code/Graphs/" *plot_name* " " *string(j)* " qubits/" * string(step) * " steps/"
        # 1 after variable names denote they're local variables in the for loop
        # And here is where the file path is defined for each iteration.
        for i = 0.00001
            title1 = "" *string(plot_name)* " plot for "* string(j) * " Qubits Jt = " * string(i)
            t_vec1, Mz_vec1 = Mz_evolve(j, step, i)
            plt=Plots.plot(t_vec1[9000:10000], Mz_vec1[9000:10000], linetype=:steppre, xlabel = "Time/ period of driving field", ylabel = " \n"*"Magnetisation / fraction of maximum value\n and orientation", legend = false)
            Plots.title!(title1)
            #Plots.savefig(figpath1*title1*".png") #Saves the plots onto my computer but requires me to make a folder

            label=title1;

            save_plot(plot_name, j, step, plt, label) # Saves the plots to github
            display(plt)

            tit1 = "FFT of "*string(plot_name)* " plot for "* string(j) * " Qubits for Jt = " * string(round(i/(10*pi), digits=2))*" π"
            Mz_vec_fft1 = fft(Mz_vec1)
            fourier=Plots.plot(abs.(Mz_vec_fft1), linetype=:steppre, xlabel="Frequecy", xlims = (0, step), ylabel="Intensity", legend = false)
            Plots.title!(tit1)
            #display(fourier)
            lab="FFT of "*string(plot_name)* " plot for "* string(j) * " Qubits for Jt = " * string(round(i/(10*pi), digits=2))*" π";

            #save_plot(plot_name, q, step, fourier, lab) # Saves the plots to github
        end
        println("Successfully finished "*string(step)*" steps\n")
    end
    println("Successfully finished "*string(j)*" qubits\n")
end

println("Successfully finished")
