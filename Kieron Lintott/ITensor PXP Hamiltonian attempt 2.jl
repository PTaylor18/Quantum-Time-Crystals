using ITensors

# This uses the ITensor system image designed to make running ITensor code
# quicker, by pre-uploading ITensor effectively
#$ julia --sysimage ~/.julia/sysimages/sys_itensors.so

function ITensors.space(
  ::SiteType"S=1/2";
  conserve_qns=false,
  conserve_sz=conserve_qns,
  conserve_szparity=false,
  qnname_sz="Sz",
  qnname_szparity="SzParity",
)
  if conserve_sz && conserve_szparity
    return [
      QN((qnname_sz, +1), (qnname_szparity, 1, 2)) => 1,
      QN((qnname_sz, -1), (qnname_szparity, 0, 2)) => 1,
    ]
  elseif conserve_sz
    return [QN(qnname_sz, +1) => 1, QN(qnname_sz, -1) => 1]
  elseif conserve_szparity
    return [QN(qnname_szparity, 1, 2) => 1, QN(qnname_szparity, 0, 2) => 1]
  end
  return 2
end

N = 30
# N = 100
S = Index(N, "Site,S=1/2")

# Arbitrary values are assigned to Up and Down state, and the Hilbert Space
# is constructed from scratch to make the Hamiltonian using operators which
# are defined below.
val(S, "Up") == 1
val(S, "Dn") == 2

ITensors.state(Up, ::SiteType"S=1/2") = [1.0, 0.0]
ITensors.state(Dn, ::SiteType"S=1/2") = [0.0, 1.0]

ITensors.state(Xp, ::SiteType"S=1/2") = [1 / sqrt(2), 1 / sqrt(2)]
ITensors.state(Xm, ::SiteType"S=1/2") = [1 / sqrt(2), -1 / sqrt(2)]

ITensors.state(Yp, ::SiteType"S=1/2") = [1 / sqrt(2), im / sqrt(2)]
ITensors.state(Ym, ::SiteType"S=1/2") = [1 / sqrt(2), -im / sqrt(2)]


# Conventional spin flip operator:
ITensors.op(::OpName"Z", ::SiteType"S=1/2") = [
  1 0
  0 -1
]

# Sz operator
ITensors.op(::OpName"Sz", ::SiteType"S=1/2") = [
  0.5 0.0
  0.0 -0.5
]

ITensors.op(::OpName"Sᶻ", t::SiteType"S=1/2") = op(OpName("Sz"), t)

ITensors.op(::OpName"S+", ::SiteType"S=1/2") = [
  0 1
  0 0
]

ITensors.op(::OpName"S⁺", t::SiteType"S=1/2") = op(OpName("S+"), t)

ITensors.op(::OpName"Splus", t::SiteType"S=1/2") = op(OpName("S+"), t)

ITensors.op(::OpName"S-", ::SiteType"S=1/2") = [
  0 0
  1 0
]

ITensors.op(::OpName"S⁻", t::SiteType"S=1/2") = op(OpName("S-"), t)

ITensors.op(::OpName"X", ::SiteType"S=1/2") = [
  0 1
  1 0
]

ITensors.op(::OpName"Sx", ::SiteType"S=1/2") = [
  0.0 0.5
  0.5 0.0
]

ITensors.op(::OpName"Sˣ", t::SiteType"S=1/2") = op(OpName("Sx"), t)

ITensors.op(::OpName"iY", ::SiteType"S=1/2") = [
  0 1
  -1 0
]

ITensors.op(::OpName"iSy", ::SiteType"S=1/2") = [
  0.0 0.5
  -0.5 0.0
]

ITensors.op(::OpName"iSʸ", t::SiteType"S=1/2") = op(OpName("iSy"), t)

ITensors.op(::OpName"Y", ::SiteType"S=1/2") = [
  0.0 -1.0im
  1.0im 0.0
]

ITensors.op(::OpName"Sy", ::SiteType"S=1/2") = [
  0.0 -0.5im
  0.5im 0.0
]

ITensors.op(::OpName"Sʸ", t::SiteType"S=1/2") = op(OpName("Sy"), t)

ITensors.op(::OpName"S2", ::SiteType"S=1/2") = [
  0.75 0.0
  0.0 0.75
]

ITensors.op(::OpName"S²", t::SiteType"S=1/2") = op(OpName("S2"), t)

# This is the all-important P operator which acts on a single site:
ITensors.op(::OpName"ProjUp", ::SiteType"S=1/2") = [
  1 0
  0 0
]

ITensors.op(::OpName"p", t::SiteType"S=1/2") = op(OpName("ProjUp"), t)

ITensors.op(::OpName"nP", ::SiteType"S=1/2") = [
  0 0
  0 1
]

ITensors.op(::OpName"np", t::SiteType"S=1/2") = op(OpName("nP"), t)

ITensors.op(::OpName"I", t::SiteType"S=1/2") = [
  1 0
  0 1
]

# With all the possible operators I'll need defined, I can start by defining
# the initial state and the Hamiltonian:
sites = siteinds("S=1/2",N; conserve_qns = false)


function hterm(sites, a, b, c)
  #=
  This function acts on a triplet of sites, and is called in later functions
  to sum over all sites to make the Hamiltonian.
  =#
  A = ITensors.op("Z", sites, a) * ITensors.op("X", sites ,b) * ITensors.op("Z", sites, c)
  B =  ITensors.op("I", sites, a) * ITensors.op("Z", sites, b) * ITensors.op("I", sites, c)
  return A + B
end

#println(hterm(sites, 1,2,3))

# Commented out section of obsolete code which may still be useful
#=
function bulksum(sites, N, periodic = true)
  bulk = Vector{ITensor}()

  for i in 1:N-2
    push!(bulk, hterm(sites, i, i+1, i+2))
  end

  if periodic == true
    push!(bulk, hterm(sites, N, 1, 2))
    push!(bulk, hterm(sites, N-1, N, 1))
  elseif periodic == false
    println("Open boundary conditions are imposed; no extra terms in H")
  else
    println("periodic must be type boolean.")
  end
  println(bulk)
  return bulk
end
=#
#=
let
  ampo = AutoMPO()
  for i = 2:N-1
    ampo += "Sz",i-1,"X",i,"Sz", i+1
    ampo += "Sz", i
  end
  H = MPO(ampo, sites)
  return H
end
=#

function Hpxp_open(sites, N)
  #=
  A function for the PXP Hamiltonian with open boundary conditions; the first
  and last site isn't considered - see notebook for how to handle these sites
  =#
  ampo = AutoMPO()
  for i = 1:N
    if i == 1
      ampo += "X",i,"Sz",i+1,"I", i+2
      ampo += "X", i, "I", i+1, "I", i+2
    elseif i == N
      ampo += "I",i-2,"Sz",i-1,"X", i
      ampo += "I", i-2, "I", i-1, "X", i
    else
      ampo += "Sz",i-1,"X",i,"Sz", i+1
      ampo += "I", i-1, "X", i, "I", i+1
    #ampo += hterm(sites, i-1, i, i+1)
    end
  end
  H = MPO(ampo, sites)
  return H, ampo
end

#println(Hpxp_open(sites, N)[1], Hpxp_open(sites, N)[2])

function Hpxp_periodic(sites, N)
  #=
  This function is the PXP Hamiltonian for periodic boundary conditions, with
  the extra measure of measuring the first and last site appropriately.
  =#
  ampo = AutoMPO()
  for i = 1:N
    if i == 1
      ampo += "Sz",N,"X",i,"Sz", i+1
      ampo += "I", N, "X", i, "I", i+1
    elseif i == N
      ampo += "Sz",N-1,"X",N,"Sz", 1
      ampo += "I", N-1, "X", N, "I", 1
    else
      ampo += "Sz",i-1,"X",i,"Sz", i+1
      ampo += "I", i-1, "X", i, "I", i+1
      #ampo += hterm(sites, i-1, i, i+1)
    end
  end
  #println(ampo)
  H = MPO(ampo, sites)
  return H, ampo
end

function Hpxp(sites, N, periodic=true)
  if periodic == true
    return Hpxp_periodic(sites, N)
  else
    return Hpxp_open(sites, N)
  end
end

state = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi0 = productMPS(sites ,state)

# These lines construct a sweeps object which is carried out 5 times
# on DMRG maxdim! sets the maximum dimensions swept for each sweep,
# and cutoff! sets the truncation error
sweeps = Sweeps(20)
maxdim!(sweeps, 10,20,100,100,200, 150, 40)
cutoff!(sweeps, 1E-10)

# This runs the DMRG algorithm included in ITensor, using psi0 as an initial
# guess for the ground state wavefunction. The optimized MPS psi and
# its eigenvalue energy are returned.
energy, psi = dmrg(Hpxp(sites, N, false)[1],psi0, sweeps)
# println(dmrg(Hpxp(sites, N, false)[1],psi0, sweeps)[1])

# With working Hamiltonians, with correct boundary conditions, the next step
# is to evolve these and see their behaviour, to ensure it does behave as
# expected compared to the PXP Hamiltonian given in the QMB Scars paper.

#=
To evolve it, you can just create the floquet operator and apply it to a state
And, I could calculate the entanglement entropy of such a state and plot that,
comparing it to the QMB Scar entanglement entropy.
=#
