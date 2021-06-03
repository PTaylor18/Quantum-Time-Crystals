using ITensors
let
  N = 100
  sites = siteinds("S=1",N)

#This tells the function Siteinds to make an array of ITensor index objects which have the properties of S=1 spins.
#This means their dimension will be 3 and they will carry the "S=1" tag,
#which will enable the next part of the code to make the appropriate operators for them

 @show sites[1]
 @show sites[2] #printing out some of these indices to verify their properties

  ampo = AutoMPO()
  for j=1:N-1
    ampo += "Sz",j,"Sz",j+1
    ampo += 1/2,"S+",j,"S-",j+1
    ampo += 1/2,"S-",j,"S+",j+1
  end
  H = MPO(ampo,sites)

#This code builds the Hamiltonian
#An AutoMPO is an object which accumulates Hamiltonian terms such as "Sz",1,"Sz",2
#so that they can be summed afterward into a matrix product operator (MPO) tensor network.
#The line of code H = MPO(ampo,sites) constructs the Hamiltonian in the MPO format, with physical indices given by the array Sites.

  psi0 = randomMPS(sites,10)

#This constructs an MPS psi0 which has the physical indices sites and a bond dimension of 10.
#It is made by a random quantum circuit that is reshaped into a MPS,
#so that it will have as generic and unbiased properties as a MPS of that size can have.
#This can help stop the DMRG calculations from getting stuck in a local minimum.


  sweeps = Sweeps(2)
  maxdim!(sweeps, 10,20,100,100,200)
  cutoff!(sweeps, 1E-10)

#Constructs a Sweeps objects which initialized to define 5 sweeps of DMRG.
#The call to macdim! Sets the maximum dimension allowed for each sweep
#The call to cutoff! Sets the truncation error goal each sweep

  energy, psi = dmrg(H,psi0, sweeps)

#The final call runs the DMRG algorithm included in ITensor, using psi0 as an initial guess for the ground state wavefunction.
#The optimized MPS psi and its eigenvalue energy are returned.

  for j=1:length(psi)
    orthogonalize!(psi,j)

# This shifts the orthogonality centre of the MPS to site number j.


    s = siteind(psi,j)

#This retrieves the site, or physical index of the jth MPS tensor.

    val = scalar(psi[j]*op(s,"Sz")*dag(prime(psi[j],s)))

#This performs the computation of the expected value of the operator "Sz".

    #println("$j $val")
  end

  b=10

  orthogonalize!(psi, b)

#	This shifts the orthogonality centre to site b of the MPS.

  U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))

#	The svd routine says to treat the link Index connecting the b'th MPS tensor psi[b]
#and the b'th physical index as "row" indicates for the purposes of the SVD

  SvN = 0.0
  for n=1:dim(S, 1)
    p = S[n,n]^2
    SvN -= p * log(p)
end
  print("SvN =",SvN)
  return
end
