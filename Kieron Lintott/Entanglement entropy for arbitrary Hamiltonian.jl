using ITensors
using Plots

N = 56
#N = 100
S = Index(N, "Site,S=1/2")
# With all the possible operators I'll need defined, I can start by defining
# the initial state and the Hamiltonian:
sites = siteinds("S=1/2",N; conserve_qns = false)

# Arbitrary values are assigned to Up and Down state, and the Hilbert Space
# is constructed from scratch to make the Hamiltonian using operators which
# are defined below.
val(S, "Up") == 1 # ground state
val(S, "Dn") == -1 # "excited" state

ITensors.state(Up, ::SiteType"S=1/2") = [1.0, 0.0]
ITensors.state(Dn, ::SiteType"S=1/2") = [0.0, 1.0]

ITensors.state(Xp, ::SiteType"S=1/2") = [1 / sqrt(2), 1 / sqrt(2)]
ITensors.state(Xm, ::SiteType"S=1/2") = [1 / sqrt(2), -1 / sqrt(2)]

ITensors.state(Yp, ::SiteType"S=1/2") = [1 / sqrt(2), im / sqrt(2)]
ITensors.state(Ym, ::SiteType"S=1/2") = [1 / sqrt(2), -im / sqrt(2)]

# Sz operator
ITensors.op(::OpName"Sz", ::SiteType"S=1/2") = [
  0.5 0.0
  0.0 -0.5
]
ITensors.op(::OpName"Z", ::SiteType"S=1/2") = [
  1.0 0.0
  0.0 -1.0
]
ITensors.op(::OpName"X", ::SiteType"S=1/2") = [
  0 1
  1 0
]
ITensors.op(::OpName"I", t::SiteType"S=1/2") = [
  1 0
  0 1
]
ITensors.op(::OpName"Is", t::SiteType"S=1/2") = [
  0.5 0.0
  0.0 0.5
]
ITensors.op(::OpName"Null", t::SiteType"S=1/2") = [
  0 0
  0 0
]

function H_N(site, N, θ)
  #=
  This function introduces the periodic driving seen in the QMB - θ = 0 is the
  default so the code runs as it did but can set it to a non-zero term to
  observe the behaviour and check the stability of time periodic behaviour
  =#
  gate = op("Is", site) - op("Sz", site)
  gates = exp(-1.0im * θ * gate)
  return gates
end

function Hpxp_open(sites, N, τ)
  #=
  A function for the PXP Hamiltonian with open boundary conditions; the first
  and last site isn't considered - see notebook for how to handle these sites
  =#
  gates = ITensor[]
  max = N-2
  for i = 1:max
    ampo = op("X",sites[i])*op("S",sites[i+1])*op("I", sites[i+2])
    ampo += op("X",sites[i])*op("I",sites[i+1])*op("I",sites[i+2])
    ampo += op("I", sites[i])*op("X", sites[i+1])*op("Z", sites[i+2])
    ampo += op("Z", sites[i])*op("X", sites[i+1])*op("I", sites[i+2])
    gatei = exp(-1.0im * τ/2 * ampo)
    push!(gates, gatei)
  end
  return gates
end

#println(Hpxp_open(sites, N)[1], Hpxp_open(sites, N)[2])

function Hpxp_periodic(sites, N, τ, θ)
  #=
  This function is the PXP Hamiltonian for periodic boundary conditions, with
  the extra measure of measuring the first and last site appropriately.
  =#
  gates = ITensor[]
  if θ == 0
    println("θ = 0, no contribution from H_N")
  else
    for i = 1:N
    #gates *= H_N(sites[i], N, θ)
    push!(gates, H_N(sites[i], N, θ))
    end
  end

  for i = 1:N
    if i == 1
      hi =  op("Z", sites[N])*op("X", sites[i])*op("Z", sites[i+1])
      hi += op("I", sites[N])*op("X", sites[i])*op("I", sites[i+1])
      hi += op("I", sites[N])*op("X", sites[i])*op("Z", sites[i+1])
      hi += op("Z", sites[N])*op("X", sites[i])*op("I", sites[i+1])
      gatei = exp(-1.0im * τ/2 * 0.5 * hi)
      #gates *= gatei
      push!(gates, gatei)
    elseif i == N
      hi =  op("Z", sites[N-1])*op("X", sites[N])*op("Z", sites[1])
      hi += op("I", sites[N-1])*op("X", sites[N])*op("I", sites[1])
      hi += op("I", sites[N-1])*op("X", sites[N])*op("Z", sites[1])
      hi += op("Z", sites[N-1])*op("X", sites[N])*op("I", sites[1])
      gatei = exp(-1.0im * τ/2 * 0.5 * hi)
      #gates *= gatei
      push!(gates, gatei)
    else
      hi =  op("Z",sites[i-1])*op("X",sites[i])*op("Z", sites[i+1])
      hi += op("I", sites[i-1])*op("X", sites[i])*op("I", sites[i+1])
      hi += op("I", sites[i-1])*op("X", sites[i])*op("Z", sites[i+1])
      hi += op("Z", sites[i-1])*op("X", sites[i])*op("I", sites[i+1])
      gatei = exp(-1.0im * τ/2 * 0.5 * hi)
      #gates *= gatei
      push!(gates, gatei)
      #ampo += hterm(sites, i-1, i, i+1)
    end
  end
  #append!(gates,reverse(gates))
  return gates
end

function Hpxp(sites, N, τ, θ, periodic=true)
  if periodic == true
    return Hpxp_periodic(sites, N, τ, θ)
  else
    return Hpxp_open(sites, N, τ)
  end
end

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

function operator(N, sites, n)
  operator = (2/N)*(op("Is", sites[n]) - op("Sz", sites[n]))
  return operator
end

function ℒ_n(psi, N, sites)
  #=
  This code aims to calculate the density imbalance of the nth site.
  =#
  #maxima = round(Int, N/2)
  DI = 0
  for i in 1:N
    orthogonalize!(psi,i)
    psidag = dag(prime(psi, "Site"))
    if isodd(i)==true
      DI += scalar(psidag[i] * operator(N, sites, i)* psi[i])
    else
      DI -= scalar(psidag[i] * operator(N, sites, i)* psi[i])
    end
  end
  #display(DI)
  return real(DI)
end

function multi_site_entropy(psi, N, sites)
  EE = 0
  Nmax = N-round(Int, N/4)
  for i in 2:Nmax
    EE += entanglement_entropy(psi, i)
  end
  #display(EE)
  return EE
end
#=
# Defines the initial state over the sites defined, which is now done in the
# function so this code is no longer required
state = [isodd(n) ? "Up" : "Up" for n=1:N]
psi0 = productMPS(sites ,state)
multi_site_entropy(psi0, N, sites)
=#
#=
# Test of ℒ_n function on a Neel State.
state = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi0 = MPS(sites ,state)

ℒ_n(psi0, N, sites)
=#
function plot_spin_and_entropy(N, ttotal, n, θ=0, periodic=true)
  # N is the number of sites, ttotal is the time the system runs over
  # and n is the site being evaluated.

  cutoff = 1E-4
  τ = 0.25

  # computes the number of steps to do
  Nsteps = Int(ttotal/τ)

  # The progress of the code is also defined
  #p = Progress(Nsteps, 1, "Computing pass ...", 50)

  # Make an array of 'site' indices
  s = siteinds("S=1/2",N;conserve_qns=false)

  # Make PXP gates, from the Hamiltonian function given above
  # This could be made into an argument of the plot function quite easily...
  gates = Hpxp(s, N, τ, θ, periodic)

  # Include gates in reverse order too
  # (N,N-1),(N-1,N-2),...
  #append!(gates,reverse(gates))

  # Initialize psi to be a product state (alternating up and down)
  psi = productMPS(s, n -> isodd(n) ? "Up" : "Dn")
  #psi = productMPS(s, n -> isodd(n) ? "Up" : "Up")

  # Compute and print initial <Sz> value and entanglement entropy of the site
  t = 0.0
  #Sz = measure_Sz(psi,n)
  #sitespin = [Sz]
  #See = [entanglement_entropy(psi, n)]
  DI = [ℒ_n(psi, N, s)]
  SEE = [multi_site_entropy(psi, N, s)]

  # This evolves for later times, where the state psi is evolved via
  # applications of gates. And Sz and e_e are evaluated for this state.
  for step=1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    #t += τ
    #Szt = measure_Sz(psi,n)
    #Seet = entanglement_entropy(psi, n)
    DIt = ℒ_n(psi, N, s)
    SEEt = multi_site_entropy(psi, N, s)
    #push!(See, Seet)
    #push!(sitespin, Szt)
    push!(DI, DIt)
    push!(SEE, SEEt)
    #println("$t $Sz")
  end

  #println(sitespin)
  #println( length(collect(0:tau:ttotal)) == length(sitespin))

  # This makes the two graphs, which are combined and outputted in a 2 graph
  # plot, which can be saved if desired
  #graphtitle = " \nlocal spin of site "*string(n)*" of "*string(N)*" sites vs time"
  #graph = Plots.plot(collect(0:τ:ttotal), sitespin, xlabel = "Time",
  #         ylabel = " \n\n"*"local spin value", title = graphtitle, size=(650, 440))
  graphtitle = " \nDensity imbalance of "*string(N)*" site system vs time"
  graph = Plots.plot(collect(0:τ:ttotal)/(1.51*π), DI, xlabel = "Time",
           ylabel = " \n \n "*"Density imbalance", title = graphtitle, size=(650, 440))
  #graph2title = "Local Entanglement Entropy of "*string(n)*" site vs time\n for a system of "*string(N)*" sites"
  #graph2 = Plots.plot(collect(0:τ:ttotal), See, xlabel = "Time",
  #         ylabel =  "\n\n"*"local entanglement entropy", title = graph2title, size=(650, 540))
  graph2title = "Entanglement Entropy of first "*string(N/2)*" sites vs time\n for a system of "*string(N)*" sites"
  graph2 = Plots.plot(collect(0:τ:ttotal)/(1.51*π), SEE, xlabel = "Time",
           ylabel =  "\n \n "*"Entanglement entropy", title = graph2title, size=(650, 540))
  twofigplot = Plots.plot(graph, graph2, layout=(2,1), show = true, legend=false, size=(700, 1020))
  #twofigtitle = "local spin and entanglement entropy of "*string(n)*" site of "*string(N)*" sites"
  twofigtitle2 = "density imbalance and entanglement entropy of first "*string(3*N/4)*" sites of "*string(N)*" sites"
  display(twofigplot)
  #savefig(twofigplot, twofigtitle*".png")

end

#0.9*π
plot_spin_and_entropy(N, 45, N-2, π , true)
