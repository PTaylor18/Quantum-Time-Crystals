using ITensors
using Tensors

import ITensors: op

N = 16
#N = 100
#S = Index(N, "Site,S=1/2")
# With all the possible operators I'll need defined, I can start by defining
# the initial state and the Hamiltonian:
sites = siteinds("S=1/2",N; conserve_qns = false)
t = 0.1;

# Sz operator, Z operator and other operators in the Hamiltonian are defined
# here, as well as the null matrix.
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
ITensors.op(::OpName"Y", ::SiteType"S=1/2") = [
  0.0 -1.0im
  1.0im 0.0
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

X = [0 1; 1 0]
I = [1 0; 0 1]

function H_N(site, N, θ)
  #=
  This function introduces the periodic driving seen in the QMB - θ = 0 is the
  default so the code runs as it did but can set it to a non-zero term to
  observe the behaviour and check the stability of time periodic behaviour
  =#
  gate = op("Is", site) - op("Sz", site)
  #gates = exp(-1.0im * θ * gate)
  return gate
end

# The sample code looked at open boundary - no (N,1) gate so I'm starting
# with the PXP hamiltonian for open boundary conditions
function Hpxp_open(sites, N, τ=0.1)
  #=
  A function for the PXP Hamiltonian with open boundary conditions; the first
  and last site isn't considered - see notebook for how to handle these sites
  =#
  gates = ITensor[]
  max = N-2
  for i = 1:max
    ampo = op("I",sites[i])*op("X",sites[i+1])*op("I", sites[i+2])
    ampo += op("I",sites[i])*op("X",sites[i+1])*op("Z",sites[i+2])
    ampo += op("Z", sites[i])*op("X", sites[i+1])*op("Z", sites[i+2])
    ampo += op("Z", sites[i])*op("X", sites[i+1])*op("I", sites[i+2])
    gatei = exp(-1.0im * τ/2 * ampo)
    # gatei = PXP_open_op(i)
    push!(gates, gatei)
  end
  display(gates)
  return gates
end

function Hpxp_open_DMRG(sites, N, τ=t)
  gate = OpSum()
  max = N-2
  for i in 1:max
    gate += "I",i,"X",i+1,"I", i+2
    gate += "I",i,"X",i+1,"Z", i+2
    gate += "Z",i,"X",i+1,"I", i+2
    gate += "Z",i,"X",i+1,"Z", i+2
  end
  return MPO(gate,sites)
end

function op(::OpName"openPXPop",::SiteType"S=1/2", s1::Index, s2::Index, s3::Index ; τ)
  ampo = op("I", s1)*op("X", s2)*op("I", s3) +
         op("I", s1)*op("X", s2)*op("Z", s3) +
         op("Z", s1)*op("X", s2)*op("Z", s3) +
         op("Z", s1)*op("X", s2)*op("I", s3)
  return exp(τ * ampo)
end
#=
function excitation_op(N, site)
  a = site - 1
  b = site + 1
  P =  [op("I", sites[1])]
  for i = 2:a
    append!(P, [op("I", sites[i])])
  end
  append!(P, [op("X", sites[site])])
  # And then apply identity operators on the rest of the sites.
  for i = b:(N-2)
    append!(P, [op("I", sites[i])])
  end
  #display(P)
  #println(typeof(P)) # An ITensor object
  return P
end
=#
function excite(N, site)
    ampo = AutoMPO()
    ampo += "S+",site
    return ampo
end
#excite(N, 4)

function deexcite(N, site)
    ampo = AutoMPO()
    ampo += "S-",site
    return ampo
end
#deexcite(N, 4)

H = Hpxp_open_DMRG(sites, N)
state = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi = productMPS(sites,state)

sweeps = Sweeps(7)
maxdim!(sweeps, 30,60,150,140,100)
cutoff!(sweeps, 1E-8)

energy, psi0 = dmrg(H,psi, sweeps)

function CdagexpCexp(ψ, s, t, N, site)
  #=
  This function creates the operator described and calculates the expectation
  value of the operator on the ground state, in order to get the Green's
  function of the system.
  ψ is the wave function, s is the site index vector, t is the time (given as a
  real quantity), N is the number of sites and site is the location of the
  excitation in the chain.
  =#
  #A = MPO(deexcite(N, site), s)
  A = op("S-", s[site])
  #A = dexcite(N, site)
  B = ops([("openPXPop", (n, n + 1, n + 2), (τ=-t * im / 2,)) for n in 1:(N - 2)], s)
  Bdag = ops([("openPXPop", (n, n + 1, n + 2), (τ=t * im / 2,)) for n in 1:(N - 2)], s)
  #Adag = MPO(excite(N, site), s)
  Adag = op("S+", s[site])

  Nsteps = 5/t
  #ρ0 = MPO(ψ)
  ρ0 = ψ
  ρ1 = ρ0
  τ = t
  for step in 1:Nsteps
    ρ1 = apply(B, ρ1)
    τ += t
  end
  #println("Got to here")
  orthogonalize!(ρ1, site)
  newA = A * ρ1[site]
  noprime!(newA)
  ρ1[site] = newA
  #println("Got to here too")
  for step in 1:Nsteps
    ρ1 = apply(Bdag, ρ1)
    #τ += t
  end
  #println("Got to here as well")
  orthogonalize!(ρ1, site)
  newAdag = Adag * ρ1[site]
  noprime!(newAdag)
  ρ1[site] = newAdag
  #println("Almost done it")
  ρ = outer(ρ1, ρ1)
  #println("\n", inner(ψ, ρ, ψ))
  return inner(ψ, ρ, ψ)
end

function expCexpCdag(ψ, s, t, N, site, totalt=5)
  #=
  This function creates the operator described and calculates the expectation
  value of the operator on the ground state, in order to get the Green's function
  of the system.
  =#
  Bdag = op("S+", s[site])
  #B = excite(N, site)
  Adag = ops([("openPXPop", (n, n + 1, n + 2), (τ=t * im / 2,)) for n in 1:(N - 2)], s)
  A = ops([("openPXPop", (n, n + 1, n + 2), (τ=-t * im / 2,)) for n in 1:(N - 2)], s)
  B = op("S-", s[site])

  Nsteps = totalt/t
  #ρ0 = MPO(ψ)
  ρ0 = ψ
  ρ1 = ρ0
  τ = t
  orthogonalize!(ρ1, site)
  newBdag = Bdag * ρ1[site]
  noprime!(newBdag)
  ρ1[site] = newBdag
  #println("step 1")
  for step in 1:Nsteps
    ρ1 = apply(A, ρ1)
    τ += t
  end
  #println("step 2")
  orthogonalize!(ρ1, site)
  newB = B * ρ1[site]
  noprime!(newB)
  ρ1[site] = newB
  #println("step 3")
  for step in 1:Nsteps
    ρ1 = apply(Adag, ρ1)
    τ += t
  end
  #println("step 4")
  ρ = outer(ρ1, ρ1)
  #println("\n", inner(ψ, ρ, ψ))
  return inner(ψ, ρ, ψ)
end


function totalGreensfunction(ψ, s, t, N, site)
  #=
  See equation referenced in my lab book. gates1 and gates2 are expCexpCdag
  and CdagexpCexp
  =#
  G = CdagexpCexp(ψ, s, t, N, site) - im * expCexpCdag(ψ, s, t, N, site)
  #println(G)
  return G
end
#=
# Make gates (1,2,3),(2,3,4),(3,4,5),...
gates = ops([("openPXPop", (n, n + 1, n + 2), (τ=-t * im / 2,)) for n in 1:(N - 2)], sites)

# This is the time evolution operator with the exponent of the Hamiltonian,
# we also want the excitation of a specific site leaving all others the same

# Include gates in reverse order too
# (N,N-1),(N-1,N-2),...
append!(gates, reverse(gates))
#display(gates)

psi = psi0
for step in 1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    τ += t
    # measure an observable
end

rho0 = outer(psi0, psi0)
rho = rho0
for step in 1:Nsteps
  rho = apply(gates, rho; cutoff=cutoff, apply_dag=true)
  τ += t
end
@show inner(psi, rho, psi)
@show inner(psi, psi)
@show tr(rho)


=#
#CdagexpCexp(psi0, sites, t, N, 3)
#totalGreensfunction(psi0, sites, t, N, 3)

Garr = [1.00 + 0.00*im]
for tl in 0.01:0.01:0.1 # tl is short for time local (to this loop in the code)
  G = totalGreensfunction(psi0, sites, tl, N, 3)
  push!(Garr, G)
end
println([0:0.01:0.1], Garr)

println("Done!")
