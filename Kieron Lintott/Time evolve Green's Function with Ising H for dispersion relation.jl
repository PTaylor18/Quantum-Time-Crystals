using ITensors
using Plots
using FFTW

import ITensors: op
#= #Defining an S=3/2 space.
ITensors.space(::SiteType"S=3/2") = 4

function ITensors.op!(Op::ITensor,
                      ::OpName"Sz",
                      ::SiteType"S=3/2",
                      s::Index)
  Op[s'=>1,s=>1] = +3/2
  Op[s'=>2,s=>2] = +1/2
  Op[s'=>3,s=>3] = -1/2
  Op[s'=>4,s=>4] = -3/2
end

function ITensors.op!(Op::ITensor,
                      ::OpName"S+",
                      ::SiteType"S=3/2",
                      s::Index)
  Op[s'=>1,s=>2] = sqrt(3)
  Op[s'=>2,s=>3] = 2
  Op[s'=>3,s=>4] = sqrt(3)
end

function ITensors.op!(Op::ITensor,
                      ::OpName"S-",
                      ::SiteType"S=3/2",
                      s::Index)
  Op[s'=>2,s=>1] = sqrt(3)
  Op[s'=>3,s=>2] = 2
  Op[s'=>4,s=>3] = sqrt(3)
end

ITensors.state(::StateName"Up", ::SiteType"S=3/2") = [1.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"Dn", ::SiteType"S=3/2") = [0.0, 0.0, 0.0, 1.0]
=# # End of the definition
N = 20
# With all the possible operators I'll need defined, I can start by defining
# the initial state and the Hamiltonian:
sites = siteinds("S=1/2",N)
t = 0.05

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

# This code relates to the PXP Hamiltonian so can be ignored.
#=
# The sample code looked at open boundary - no (N,1) gate so I'm starting
# with the PXP hamiltonian for open boundary conditions.

function Hpxp_open(sites, N, τ=t)
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
=#
function IsingH_DMRG(sites, N, S)
  # sites is the 1D array of spin sites, N is the number of sites and
  # S is the interaction strength between nearest neighbours.
  ampo = OpSum()
  max = N - 1
  for j=1: max
    ampo += S,"Sz",j,"Sz",j+1
    ampo += 1/2,"Sz",j
  end
  ampo += 1/2,"Sz",N
  H = MPO(ampo,sites)
  return H
end

function op(::OpName"openIsingH", ::SiteType"S=1/2", s1::Index, s2::Index; τ, S)
  ampo = S * op("Sz", s1)*op("Sz", s2)
       + 0.5 * op("Sz", s1)
  return exp(ampo * τ)
end

function op(::OpName"Ck", ::SiteType"S=3/2", s1::Index; C, N, nk)
  # S+ should have +ve nk and S- should have -ve nk, as S- corresponds to the
  # creation of a quasiparticle and the creation operator has the -ve sign.
  coeff = exp(im * (2*pi)^2 * nk * 1/N) #* N^(-0.5)
  ampoj = coeff * op(C, s1)
  return ampoj
end

function Cop(ρ1, site, C)
  #=
  This operator takes the operator C, which may be the excitation or the
  deexcitation operator and applies it to the state ρ1 acting on the site
  number "site" of the N site system.
  In principle, this could be applied to different sites to simulate multiple
  excitations.
  =#
  orthogonalize!(ρ1, site)
  newC = C * ρ1[site]
  noprime!(newC)
  ρ1[site] = newC
  return ρ1
end

function Ckop(s, C, N, nk)
  #=
  This function is Ck, for k not equal to zero.
  ρ1 and C are the state and form of the single site operator,
  N is the total number of sites, nk is the wavenumber and

  Ckop is called in the place of Cop, so the operator acts on a given
  k-state and the Green's function will be a function of the wavevector
  of the excitation in the system.
  To handle the complex conjugate for the other operator, the sign is
  changed by the value of nk passed.
  =#
  ampo = ITensor[]
  #min = -N/2; max = abs(min)

  for j in 1:N
    sj = s[j]
    coeff = exp(im * (2*pi)^2 * nk * j/N) #* N^0.5
    ampoj = coeff * op(C, sj)
    #println(typeof(ampoj))
    push!(ampo, ampoj)
  end
  #println(typeof(ampo))
  return ampo
  #return MPO(ampo, s)
end

function Ck_MPO(sites, C, N, nk, S=-1)
  # sites is the 1D array of spin sites, N is the number of sites and
  # S is the interaction strength between nearest neighbours.
  ampo = OpSum()
  for j = 1: N
    coeff = S * exp(im * (2*pi) * nk * j/N)
    ampo += coeff, C, j
  end
  O = MPO(ampo,sites)
  return O
end

function Getψ(sites, N, J)
  # Computing the ground state wavefunction of the PXP Hamiltonian
  # using the DMRG algorithm
  H = IsingH_DMRG(sites, N, J)
  state = [isodd(n) ? "Up" : "Up" for n=1:N]
  psi = productMPS(sites,state)

  sweeps = Sweeps(5)
  maxdim!(sweeps, 30,60,150,100,40)
  cutoff!(sweeps, 1E-8)

  println("Energies for J = "*string(J)*" :\n ")
  energy, psi0 = dmrg(H, psi, sweeps)
  println("\n")

  return psi0
end

function CdagexpCexp(ψ, s, t, N, totalt, S, nk)
  # The number of sites in the 1D chain
  #=
  This function creates the operator described and calculates the expectation
  value of the operator on the ground state, in order to get the Green's
  function of the system which compromises of this value and another similar
  #
  value due to the next function defined, expCexpCdag().
  ψ is the wave function, s is the site index vector, t is the time (given as a
  # Sz operator, Z operator and other operators in the Hamiltonian are defined
  real quantity), N is the number of sites and site is the location of the
  excitation in the chain.
  S is nearest neighbour interaction strength
  =#

  #sid = siteind(ψ, site)
  A = ops([("Ck", n, ("S+", N, nk, )) for n in 1:N], s)
  B = ops([("openIsingH", (n, n + 1), (τ=-t * im, S, )) for n in 1:(N - 1)], s)
  Bdag = ops([("openIsingH", (n, n + 1), (τ=t * im, S, )) for n in 1:(N - 1)], s)
  Adag = ops([("Ck", n, ("S-", N, -nk, )) for n in 1:N], s)

  Nsteps = totalt/t
  #ρ0 = MPO(ψ)
  #ρ0 = ψ
  ρ1 = copy(ψ)
  #println(τ)
  τl = 0.0
  for step in 1:Nsteps
    ρ1 = apply(B, ρ1; cutoff=1E-12)
    #τl += t/2
    #println(τl)
  end
  #println("Got to here")
  #Old code pertaining to the de-excitation operator
  ρ1 = apply(A, ρ1; cutoff=1E-12)

  #ρ1 = deexcite(s, site) *   ρ1
  #println("Got to here too")
  for step in 1:Nsteps
    ρ1 = apply(Bdag, ρ1; cutoff=1E-12)
  end
  #println("Got to here as well")
  #Old code pertaining to the excitation operator
  #=orthogonalize!(ρ1, site)
  newAdag = Adag * ρ1[site]
  noprime!(newAdag)
  ρ1[site] = newAdag=#
  ρ1 = apply(Adag, ρ1; cutoff=1E-12)
  #println("Almost done it")
  #println("\n", inner(ψ, ρ1))
  return inner(ψ, ρ1)/N
end

function expCexpCdag(ψ, s, t, N, totalt, S, nk)
  #=
  This function creates the operator described and calculates the expectation
  value of the operator on the ground state, in order to get the Green's function
  of the system.
  =#
  #sid = siteind(ψ, site)
  Bdag = ops([("Ck", n, ("S-", N, -nk, )) for n in 1:N], s)
  B = ops([("Ck", n, ("S+", N, nk, )) for n in 1:N], s)
  Adag = ops([("openIsingH", (n, n + 1), (τ=t * im, S,)) for n in 1:(N - 1)], s)
  A = ops([("openIsingH", (n, n + 1), (τ=-t * im, S,)) for n in 1:(N - 1)], s)
  #B = op(sid, "S-")

  Nsteps = totalt/t
  #ρ0 = MPO(ψ)
  #ρ0 = ψ
  ρ1 = copy(ψ)
  #=
  orthogonalize!(ρ1, site)
  newBdag = Bdag * ρ1[site]
  noprime!(newBdag)
  ρ1[site] = newBdag=#
  #ρ1 = Cop(ρ1, site, Bdag)
  ρ1 = apply(Bdag, ρ1; cutoff=1E-12)
  #println("step 1")
  for step in 1:Nsteps
    ρ1 = apply(A, ρ1; cutoff=1E-12)
    #τ += t
  end
  #println("step 2")
  #ρ1 = Cop(ρ1, site, B)
  ρ1 = apply(B, ρ1; cutoff=1E-12)
  #println("step 3")
  for step in 1:Nsteps
    ρ1 = apply(Adag, ρ1; cutoff=1E-12)
    #τ += t
  end
  #println("step 4")
  #ρ = outer(ρ1, ρ1)
  #println("\n", inner(ψ, ρ1))
  return inner(ψ, ρ1)/N
end


function totalGreensfunction(ψ, s, t, N, totalt, S, nk)
  #=
  See equation referenced in my lab book. gates1 and gates2 are expCexpCdag
  and CdagexpCexp
  =#
  Val = CdagexpCexp(ψ, s, t, N, totalt, S, nk) + expCexpCdag(ψ, s, t, N, totalt, S, nk)
  #println(G)
  G = -im * Val
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

#CdagexpCexp(psi0, sites, 0.05, N, 3, 0.4)
#totalGreensfunction(psi0, sites, 0.001, N, 3, 1.3)

function SpectralFunction(X, A, N, St, nk)
  #=
  A is a 1D array with values of frequency in X and the Fourier Transformed
  Green's Function in A, so Y is the Spectral Function.
  =#
  Y = - (1/pi) * A
  display(plot(X, abs.(Y), xlabel = "frequency\n ", ylabel = " \n"*"Spectral Function",
                title=" \nSpectral Function vs frequency for "*string(N)*" sites,\n wavevector = "*string(nk)*" and J = "*string(St)*"\n ",
                legend = false))
end

function FFT_plan(tmax, t)
  # Making the FFT plan. tmax, t are variables defined in another function
  len = 1 + tmax / t
  X = rand(ComplexF64, round(Int, len))
  P = plan_fft(X)
  return P
end
# And then I think P can be used in place.

function PlotG(sites, t, N, tmax, S, nk)

  psi0 = Getψ(sites, N, S)

  Garrc = [0.00 - im*1.00];
  for tl in t:t:tmax # tl is short for time local (to this loop in the code)
    G = totalGreensfunction(psi0, sites, t, N, tl, S, nk)
    #push!(Garr, real(G))
    push!(Garrc, G)
    #println(tl)
  end
  #println([0:0.01:1], Garr)

  orig_t = 0:t:tmax
  sr = 1.0/t
  ωx = fftfreq(length(orig_t), sr) |> fftshift
  #println(ωx)

  P = FFT_plan(tmax, t)
  Yp = P * Garrc
  # Yp could then be used in place of Y

  #Y = fft(Garrc)
  #iY = imag.(Y)
  Yshift = fftshift(Yp)
  #=
  println("frequency array length is: ", length(ωshift),"\n")
  println("Fourier component array length is: ",length(Yshift))
  =#
  header = " \nGreen's Function vs time for "*string(N)*" sites\n in the system\n"
  display(plot(orig_t, imag.(Garrc), title=header, legend = false))
  #println("Green's Function plot made.")
  #=
  FTplot = plot(abs.(Yshift), xlabel = "frequency", ylabel = " \n"*"Fourier components",
                title=" \nAbsolute value of the Fourier Transform\n of Green's Function vs frequency for "*string(N)*" sites \nwith site "*string(site)*" excited\n \n", legend = false)
  display(FTplot)
  =#
  Ysf = imag.(Yshift) # This is the imaginary part of the FT of the Retarded
                              # Green's Function needed to compute the Spectral Function.
  #Ysf = imag(fft(Garrc))
  #println(Ysf)
  SpectralFunction(ωx, Ysf, N, S, nk)
end

# Some more sample code about applying the excitation operators in the most
# effective way; I think this approach is the best that can be done...
#=
j = 4 # called site in my code
psi = copy(psi0)
s = siteind(psi,j)
newpsi= op(s,"S+") * psi[j]
noprime!(newpsi)
psi[j]= newpsi
println("First value is ", inner(psi, psi0))
=#
#PlotG(sites, 0.05, N, 15, -1, 1)


# A test of the quasiparticle operators, as it didn't quite work the first time
nk = 0
psi = Getψ(sites, N, -1)
#psi1 = copy(psi)
O = Ck_MPO(sites, "S+", N, nk)
#Os = op("S-", sites[2])
#Osdag = op("S+", sites[2])
#O1 = ops([("Ck", (n), ("S+", N, nk, ))  for n in 1:(N)], sites)
Odag = Ck_MPO(sites, "S-", N, -nk)
#O1dag = ops([("Ck", (n), ("S+", N, nk, ))  for n in 1:(N)], sites)
res = applyMPO(O, psi)
res1 = applyMPO(Odag, res)
display(inner(psi1, res1))

#ops([("openIsingH", (n, n + 1), (τ=t * im, S,)) for n in 1:(N - 1)], s)

println("Done!\n ")
#=
t = 0/05; tmax = 5
orig_t = 0:t:tmax
sr = 1.0/t
ωx = fftfreq(length(orig_t), sr) |> fftshift
println(ωx)

tmat = 0:t:tmax
freqs = fftfreq(length(tmat), 1.0/Ts) |> fftshift
println(freqs[1:10])
=#
