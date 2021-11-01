using ITensors

N = 14
#N = 100
#S = Index(N, "Site,S=1/2")
# With all the possible operators I'll need defined, I can start by defining
# the initial state and the Hamiltonian:
sites = siteinds("S=1/2",N; conserve_qns = false)

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
Y = [0 -im; im 0]
Z = [1 0; 0 -1]
# Define the Hamiltonian you want to construct.
# In my case this is the PXP hamiltonian which requires
# a rotation term, and operators acting on 3 sites which
# differ by their periodic boundary condition, hence why
# this is handed by Hpxp which calls one of the functions
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

# H_n isn't called in this function as we don't consider this term in the open
# boundary conditions just yet.
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

# H_n is called in this function
function Hpxp_periodic(sites, N, τ, θ)
  #=
  This function is the PXP Hamiltonian for periodic boundary conditions, with
  the extra measure of measuring the first and last site appropriately.
  =#
  gates = AutoMPO()
  if θ == 0
    println("θ = 0, no contribution from H_N")
  else
    for i = 1:N
    #gate *= H_N(sites[i], N, θ)
    #gates += H_N(sites[i], N, θ)
    gates += "Is",i
    gates += -1.0,"Sz",i
    end
  end

  for i = 1:N
    if i == 1
      #=
      gates += op("Z", sites[N])*op("X", sites[i])*op("Z", sites[i+1])
      gates += op("I", sites[N])*op("X", sites[i])*op("I", sites[i+1])
      gates += op("I", sites[N])*op("X", sites[i])*op("Z", sites[i+1])
      gates += op("Z", sites[N])*op("X", sites[i])*op("I", sites[i+1])
      gatei = exp(-1.0im * τ/2 * 0.5 * hi)
      gates *= gatei
      push!(gates, gatei)
      =#
      gates += "Z",N,"X",i,"Z",i+1
      gates += "I",N,"X",i,"I",i+1
      gates += "I",N,"X",i,"Z",i+1
      gates += "Z",N,"X",i,"I",i+1
    elseif i == N
      #=
      gates += op("Z", sites[N-1])*op("X", sites[N])*op("Z", sites[1])
      gates += op("I", sites[N-1])*op("X", sites[N])*op("I", sites[1])
      gates += op("I", sites[N-1])*op("X", sites[N])*op("Z", sites[1])
      gates += op("Z", sites[N-1])*op("X", sites[N])*op("I", sites[1])
      gatei = exp(-1.0im * τ/2 * 0.5 * hi)
      gates *= gatei
      push!(gates, gatei)
      =#
      gates += "Z",N-1,"X",N,"Z",1
      gates += "I",N-1,"X",N,"I",1
      gates += "I",N-1,"X",N,"Z",1
      gates += "Z",N-1,"X",N,"I",1
    else
      #=
      gates += op("Z",sites[i-1])*op("X",sites[i])*op("Z", sites[i+1])
      gates += op("I", sites[i-1])*op("X", sites[i])*op("I", sites[i+1])
      gates += op("I", sites[i-1])*op("X", sites[i])*op("Z", sites[i+1])
      gates += op("Z", sites[i-1])*op("X", sites[i])*op("I", sites[i+1])
      gatei = exp(-1.0im * τ/2 * 0.5 * hi)
      gates *= gatei
      push!(gates, gatei)
      =#
      gates += "Z",i-1,"X",i,"Z",i+1
      gates += "I",i-1,"X",i,"I",i+1
      gates += "I",i-1,"X",i,"Z",i+1
      gates += "Z",i-1,"X",i,"I",i+1
    end
  end
  #append!(gates,reverse(gates))

  return MPO(gates, sites)
end

# This function differentiates between the boundary conditions by the boolean
# key word argument in the function.
function Hpxp(sites, N, τ, θ, periodic=true)
  if periodic == true
    return Hpxp_periodic(sites, N, τ, θ)
  else
    return Hpxp_open(sites, N, τ)
  end
end

# Use DMRG to find the ground state of a given Hamiltonian
τ = 0.2; θ = 0.01; #Hpxp(sites, N, τ, θ, true)
H = Hpxp(sites, N, τ, θ, true)
state = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi0 = productMPS(sites,state)

sweeps = Sweeps(7)
maxdim!(sweeps, 30,60,150,140,100)
cutoff!(sweeps, 1E-8)

energy, psi = dmrg(H,psi0, sweeps)

# Use TEBD to evolve the Hamiltonian in time.

# Maybe define the a and a† operators needed for the correlation
# function - will have to modify the output of the functions to get the
# exponential factor.
#= The function can be run with the other code block I commented out from the
previous code to calculate the gates which can be used for TEBD in the usual
way, unless MPO Hamiltonian can be used in TEBD.
=#

# Compute Pauli matrices tensor product.
#ITensors.op(::OpName"A", ::SiteType"S=1/2") =
#M1 = kron(Y, Z)
#display(kron(X, M1))
function ITensors.op(::OpName"PauliTP", ::SiteType"S=1/2",
   s1::Index, s2::Index)
  M1 = kron(Y, Z)
  M2 = kron(X, M1)
  return M2
end

# Here I define the operator needed to compute the correlation function

# Evaluate the correlation function of the wavefunction psi at
# each step of the TEBD code?

# Sum over all possible terms to get the Green's function.
