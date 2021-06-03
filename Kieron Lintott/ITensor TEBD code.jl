using ITensors
using Plots

let
  N = 150
  cutoff = 1E-8
  tau = 0.1
  ttotal = 5.0
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

  # Function that measures <Sz> on site n
  function measure_Sz(psi,n)
    psi = orthogonalize(psi,n)
    sn = siteind(psi,n)
    Sz = scalar(dag(prime(psi[n],"Site"))*op("Sz",sn)*psi[n])
    return real(Sz)
  end

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

  # This gives the site in the middle of the chain on N sites
  c = round(div(N,2)+2)
  #c = 23

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
    push!(sitespin, Szt)
    #println("$t $Sz")
  end

  #println(sitespin)
  #println( length(collect(0:tau:ttotal)) == length(sitespin))
  graphtitle = "local spin of site "*string(c)*" of "*string(N)*" sites over time"
  graph = Plots.plot(collect(0:tau:ttotal), sitespin, xlabel = "Time",
           ylabel = "local spin value", title = graphtitle)
  display(graph)
  savefig(graph, graphtitle*".png")

  return
end
