using ITensors
using Plots
theme(:dao)

function save_plot(name, N, plot, label)
    cd(dirname(@__FILE__));
    dir = pwd();
    println("Saving plots to the directory in $dir","\\Graphs")
    strlabel = string(label);
    mkpath(string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",N, " sites"))
    pngfilename = string(dir,"\\Graphs\\","\\",name,"\\","\\",name,"_",N, " sites\\",strlabel,".png")
    savefig(plot, pngfilename)
end

macro Name(arg)
   string(arg)
end

plot_name = @Name Local_site


let
  N = 150
  cutoff = 1E-8
  tau = 0.1
  ttotal = 10.0
  println("Successfully started \n")
  # Compute the number of steps to do
  Nsteps = Int(ttotal/tau)

  # Make an array of 'site' indices
  # Defiens array of spin 1/2 tensor indices
  s = siteinds("S=1/2",N;conserve_qns=true)

  # Make gates (1,2),(2,3),(3,4),...
  # Makes an empty array that will hold ITensors
  # that will be Trotter gates
  gates = ITensor[]
  for j=1:N-1
    s1 = s[j]
    s2 = s[j+1]
    hj =       op("Sz",s1) * op("Sz",s2) +
         1/2 * op("S+",s1) * op("S-",s2) +
         1/2 * op("S-",s1) * op("S+",s2)
    # The op function reads "S=1/2" tag on the site indices (sites j and j+1)

    Gj = exp(-1.0im * tau/2 * hj)
    # Makes the Trotter gate

    push!(gates,Gj)
  end
  # Include gates in reverse order too
  # to complete the correct trotter formula
  # (N,N-1),(N-1,N-2),...
  append!(gates,reverse(gates))

  # Function that measures <Sz> on site n
  function measure_Sz(psi,n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi,n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
    return real(Sz)
  end

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  c = div(N,2)

  # Compute and print initial <Sz> value
  t = 0.0
  Sz = measure_Sz(psi,c)
  #println("$t $Sz")

  # Do the time evolution by applying the gates
  # for Nsteps steps
  spin = [Sz]
  for step=1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    t += tau
    Sz1 = measure_Sz(psi,c)
    push!(spin, Sz1)
    #println("$t $Sz")
  end

  title1 = "Local spin of site "*string(c)*" of "*string(N)*" sites over time"
  fig = Plots.plot(0:tau:ttotal, spin, xlabel= "Time", ylabel= "Spin", title=title1,  legend = false)
  Plots.title!(title1)
  label = title1;

  display(fig)
  save_plot(plot_name, N, fig, label) # Saves the plots to github

  println("Successfully finished \n")
  return
end
