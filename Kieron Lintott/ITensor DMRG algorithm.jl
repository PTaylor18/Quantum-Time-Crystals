using ITensors
let
  #These lines make the array with S=1 spins. It's defined in this
  # way so the following code can define appropriate operators
  N = 100
  sites = siteinds("S=1",N; conserve_qns = true)

  # The addition of the kwarg = true is that every qn associated with every
  # site is conserved. There're other keywords which can conserve a subset
  # of a QN. This check shows QNs are conserved:
  # @show sites[3]

  # This builds the Hamiltonian of the system we want to describe, a 1D
  # Heisenberg model of integer spins.
  ampo = AutoMPO()
  for j=1:N-1
    ampo += "Sz",j,"Sz",j+1
    ampo += 1/2,"S+",j,"S-",j+1
    ampo += 1/2,"S-",j,"S+",j+1
  end
  # These terms above are compiled into a Matrix Product Operator
  # tensor network which sums the individual contributions over
  # all the sites in the sites array.
  H = MPO(ampo,sites)

  # This line constructs a matrix product state with the indices
  # from sites and a bond dimension of 10. A random Quantum Circuit
  # is reshaped into an MPS.
  # psi0 = randomMPS(sites,10)
  # These changes are made for the conserved QNs
  state = [isodd(n) ? "Up" : "Dn" for n=1:N]
  psi0 = productMPS(sites,state)

  # This makes an alternating spin state with neighbouring sites having
  # antiparallel spins. The second line takes the state array and
  # applies it to the sites array. Such a state is a Neel state, with a total
  # spin of zero (for an even number of sites).

  # Alternatively could make a state array with up, dn, Z0 representing 1, 0,
  # -1 spin states for S = 1 particles.

  # These lines construct a sweeps object which is carried out 5 times
  # on DMRG maxdim! sets the maximum dimensions swept for each sweep,
  # and cutoff! sets the truncation error
  sweeps = Sweeps(6)
  maxdim!(sweeps, 10,20,100,100,200, 150)
  cutoff!(sweeps, 1E-10)

  # This runs the DMRG algorithm included in ITensor, using psi0 as an initial
  # guess for the ground state wavefunction. The optimized MPS psi and
  # its eigenvalue energy are returned.
  energy, psi = dmrg(H,psi0, sweeps)

  #println(psi)

  # A common task is to measure local observables (at a given site) of a
  # wavefunction, which can be easily done by a for loop over all the sites
  # given in the MPS:

  for j=1:length(psi)
    # See comments on the website, but it effectively only uses the tensor at site
    # j as the tensors from all other sites don't contribute to the expectation
    # values at that site.
    orthogonalize!(psi,j)

    # This retrieves the site/ tensor index. Could also obtain it from the
    # siteinds array calculated earlier.
    s = siteind(psi,j)
    # This next line calculates the expectation value at the site. The call to
    # op(s,"Sz") makes the ITensor for the "Sz" operator, where here we are
    # assuming that the Index s carries a physical tag type such that "Sz" is
    # defined for this Index. (Examples could include the tag "S=1/2" or the tag
    # "Electron".) The returned operator has two indices: s and s'. Contracting
    # the operator with psi[j] contracts over the Index s. Then contracting with
    # dag(prime(psi[j],s)) contracts over s' and the two bond indices of the MPS
    # conecting to the jth MPS tensor. Finally, the call to scalar converts the
    # resulting scalar-valued ITensor into a number which is stored in val:
    val = scalar(psi[j]*op(s,"Sz")*dag(prime(psi[j],s)))

    # Prints out val of jth site, could easily store it in an array to
    # plot as a graph if suitable for the relevant purpose
    #println("$j $val")
  end

  # Calculating entanglement entropy of psi state
  b = 33
  orthogonalize!(psi, b)
  U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
  SvN = 0.0
  for n=1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
  end
  print(SvN)
  return
end

# you can take the returned MPS psi and do further calculations with it,
# such as measuring local operators or computing entanglement entropy.

# Results for original code as given on ITensor website:
# 1 energy=-138.822059494819 maxlinkdim=10 maxerr=1.41E-02 time=34.537
# After sweep 2 energy=-138.937257242465 maxlinkdim=20 maxerr=4.86E-06 time=2.164
# After sweep 3 energy=-138.940084410127 maxlinkdim=89 maxerr=9.97E-11 time=8.828
# After sweep 4 energy=-138.940086057246 maxlinkdim=99 maxerr=9.99E-11 time=17.814
# After sweep 5 energy=-138.940086056070 maxlinkdim=95 maxerr=9.99E-11 time=18.425
