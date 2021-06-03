using ITensors

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
  println("\nEntanglement entropy is: ", SvN, "\n")
  return SvN
end

function Psi(N::Integer, b)
  # This function makes an arbitrary state which we calculate
  # the entanglement entropy of in the next step.
  let
    sites = siteinds("S=1",N)

    ampo = AutoMPO()
    for j=1:N-1
      ampo += "Sz",j,"Sz",j+1
      ampo += 1/2,"S+",j,"S-",j+1
      ampo += 1/2,"S-",j,"S+",j+1
    end
    H = MPO(ampo,sites)

    psi0 = randomMPS(sites,10)

    sweeps = Sweeps(5)
    maxdim!(sweeps, 10,20,100,100,200)
    cutoff!(sweeps, 1E-10)

    energy, psi_local = dmrg(H,psi0, sweeps)

    entanglement_entropy(psi_local, b)
  end
  #println(psi)
end

Psi(100, 9)

# N.B. Get positive entangelement entropies despite the
# negative sign in the sum. I think this is because the
# system is losing its entanglement so the entropy is
# increasing to indicate disorder in the system...

# N.B. Found entangelement entropy is higher in the bulk
# than on the edges.
