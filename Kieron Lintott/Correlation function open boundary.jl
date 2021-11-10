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

function excitation_op(N, site)
  a = site - 1
  b = site + 1
  P = otimes(op("I", sites[2]), op("I", sites[1]))
  for i = 3:a
    P = otimes(op("I", sites[i]), P)
  end
  Px = otimes(op("X", sites[site]), P)
  # And then apply identity operators on the rest of the sites.
  for i = b:N
    Px = otimes(op("I", sites[i]), P)
  end
  display(Px)
  #println(typeof(Px)) # An ITensor object
  return Px
end

H = Hpxp_open_DMRG(sites, N)
state = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi0 = productMPS(sites,state)

sweeps = Sweeps(7)
maxdim!(sweeps, 30,60,150,140,100)
cutoff!(sweeps, 1E-8)

energy, psi1 = dmrg(H,psi0, sweeps)

Nsteps = 5/t

function CdagexpCexp(ψ, s, t, N, site)
  #=
  This function creates the operator described and calculates the expectation
  value of the operator on the ground state, in order to get the Green's
  function of the system.
  ψ is the wave function, s is the site index vector, t is the time (given as a
  real quantity), N is the number of sites and site is the location of the
  excitation in the chain.
  =#
  A = excitation_op(N, site)
  B = ops([("openPXPop", (n, n + 1, n + 2), (τ=-t * im / 2,)) for n in 1:(N - 2)], s)
  gates = A * B
  append!(gates, reverse(gates))
  return gates
end

function expCexpCdag(ψ, s, t)
#=
This function creates the operator described and calculates the expectation
value of the operator on the ground state, in order to get the Green's function
of the system.
=#
  B = excitation_op(N, site)
  A = ops([("openPXPop", (n, n + 1, n + 2), (τ=t * im / 2,)) for n in 1:(N - 2)], s)
  gates = A * B
  append!(gates, reverse(gates))
  return gates
end
#=
# Make gates (1,2,3),(2,3,4),(3,4,5),...
gates = ops([("openPXPop", (n, n + 1, n + 2), (τ=-t * im / 2,)) for n in 1:(N - 2)], sites)

# This is the time evolution operator with the exponent of the Hamiltonian,
# we also want the excitation of

# Include gates in reverse order too
# (N,N-1),(N-1,N-2),...
append!(gates, reverse(gates))
#display(gates)
# Note: I think you need to

psi = psi1
for step in 1:Nsteps
    psi = apply(gates, psi; cutoff=cutoff)
    τ += t
    # measure an observable
end

rho0 = outer(psi1, psi1)
rho = rho0
for step in 1:Nsteps
  rho = apply(gates, rho; cutoff=cutoff, apply_dag=true)
  τ += t
end
@show inner(psi, rho, psi)
@show inner(psi, psi)
@show tr(rho)


=#
CdagexpCexp(psi1, sites, t, N, 5)
println("Done!")
