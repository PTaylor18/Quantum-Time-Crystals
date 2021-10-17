using ITensors
using Plots
using ProgressMeter

function entanglement_entropy(psi, b)
  # This function calculates the entanglement entropy of a state
  # of N sites for the b-th site.
  #println(psi)
  orthogonalize!(psi, b)
  U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
  SvN = 0.0
  for n=1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
  end
  #println("\nEntanglement entropy is: ", SvN, "\n")
  return SvN
end

function measure_Sz(psi,n)
  # Function that measures <Sz> on site n for state psi.
  psi = orthogonalize(psi,n)
  sn = siteind(psi,n)
  Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
  return real(Sz)
end

function plot_spin_and_entropy(N, ttotal, n)
  # N is the number of sites, ttotal is the time the system runs over
  # and n is the site being evaluated.

  cutoff = 1E-5
  tau = 0.2

  # computes the number of steps to do
  Nsteps = Int(ttotal/tau)

  # The progress of the code is also defined
  #p = Progress(Nsteps, 1, "Computing pass ...", 50)

  # Make an array of 'site' indices
  s = siteinds("S=1/2",N;conserve_qns=true)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ITensor[]
  for j=1:N-1
    s1 = s[j]
    s2 = s[j+1]
    # hj is the operator which is summed over j to return the total
    # Hamiltonian, known as the 1D Heisenberg Hamiltonian.
    hj =       op("Sz",s1) * op("Sz",s2) +
         1/2 * op("S+",s1) * op("S-",s2) +
         1/2 * op("S-",s1) * op("S+",s2)
    # Trotter gate Gj is appened to the gates array.
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
  end
  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  append!(gates,reverse(gates))

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  # This gives the site in the middle of the chain on N sites
  #c = round(div(N,2)+1)
  # Compute and print initial <Sz> value and entanglement entropy of the site
  t = 0.0
  Sz = measure_Sz(psi,n)
  sitespin = [Sz]
  See = [entanglement_entropy(psi, n)]

  # This evolves for later times, where the state psi is evolved via
  # applications of gates. And Sz and e_e are evaluated for this state.
  for step=1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    t += tau
    Szt = measure_Sz(psi,n)
    Seet = entanglement_entropy(psi, n)
    push!(See, Seet)
    push!(sitespin, Szt)
    #println("$t $Sz")
  end

  #println(sitespin)
  #println( length(collect(0:tau:ttotal)) == length(sitespin))

  # This makes the two graphs, which are combined and outputted in a 2 graph
  # plot, which can be saved if desired
  graphtitle = " \nlocal spin of site "*string(n)*" of "*string(N)*" sites vs time"
  graph = Plots.plot(collect(0:tau:ttotal), sitespin, xlabel = "Time",
           ylabel = " \n\n"*"local spin value", title = graphtitle, size=(650, 440))
  graph2title = "Local Entanglement Entropy of "*string(n)*" site vs time\n for a system of "*string(N)*" sites"
  graph2 = Plots.plot(collect(0:tau:ttotal), See, xlabel = "Time",
           ylabel = "local entanglement entropy", title = graph2title, size=(650, 540))
  twofigplot = Plots.plot(graph, graph2, layout=(2,1), show = true, legend=false, size=(700, 1020))
  twofigtitle = "local spin and entanglement entropy of "*string(n)*" site of "*string(N)*" sites"
  display(twofigplot)
  #savefig(twofigplot, twofigtitle*".png")

end

N=30
plot_spin_and_entropy(N, 10, 14)

#= Note:

=#

#=
# This section of code was trial code used to write the above function
let
  N = 100
  cutoff = 1E-7
  tau = 0.1
  ttotal = 10.0
  #println("This is running")

  # Compute the number of steps to do
  Nsteps = Int(ttotal/tau)

  # Make an array of 'site' indices
  s = siteinds("S=1/2",N;conserve_qns=true)

  # Make gates (1,2),(2,3),(3,4),...
  gates = ITensor[]
  for j=1:N-1
    s1 = s[j]
    s2 = s[j+1]
    # hj is the operator which is summed over j to return the total
    # Hamiltonian.
    hj =       op("Sz",s1) * op("Sz",s2) +
         1/2 * op("S+",s1) * op("S-",s2) +
         1/2 * op("S-",s1) * op("S+",s2)
    # Trotter gate Gj is appened to the gates array.
    Gj = exp(-1.0im * tau/2 * hj)
    push!(gates,Gj)
  end
  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  append!(gates,reverse(gates))

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  # This gives the site in the middle of the chain on N sites
  #c = round(div(N,2)+1)
  c = N-1
  See = [entanglement_entropy(psi, c)]

  # Compute and print initial <Sz> value
  t = 0.0
  Sz = measure_Sz(psi,c)
  #println("$t $Sz")

  #println("This is still running")
  # Do the time evolution by applying the gates
  # for Nsteps steps
  sitespin = [Sz]
  for step=1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    t += tau
    Szt = measure_Sz(psi,c)
    Seet = entanglement_entropy(psi, c)
    push!(See, Seet)
    push!(sitespin, Szt)
    #println("$t $Sz")
  end

  #println(sitespin)
  #println( length(collect(0:tau:ttotal)) == length(sitespin))
  graphtitle = "local spin of site "*string(c)*" of "*string(N)*" sites vs time"
  graph = Plots.plot(collect(0:tau:ttotal), sitespin, xlabel = "Time",
           ylabel = " \n\n"*"local spin value", title = graphtitle, size=(650, 460))
  graph2title = "Local Entanglement Entropy of "*string(c)*" site vs time\n for a system of "*string(N)*" sites"
  graph2 = Plots.plot(collect(0:tau:ttotal), See, xlabel = "Time",
           ylabel = "local entanglement entropy", title = graph2title, size=(650, 540))
  twofigplot = Plots.plot(graph, graph2, layout=(2,1), legend=false, size=(700, 1020))
  display(twofigplot)
  twofigtitle = "local spin and entanglement entropy of "*string(c)*" site of "*string(N)*" sites"
  #savefig(twofigplot, twofigtitle*".png")

  return
end
=#
