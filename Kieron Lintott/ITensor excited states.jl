using ITensors

let
  # Size of the system being simulated
  N = 20

  sites = siteinds("S=1/2",N)

  h = 4.0

  weight = 10*h # use a large weight
                # since gap is expected to be large

  #
  # Use the AutoMPO feature to create the
  # transverse field Ising model
  #
  # Factors of 4 and 2 are to rescale
  # spin operators into Pauli matrices
  #
  ampo = AutoMPO()
  for j=1:N-1
    ampo += -4,"Sz",j,"Sz",j+1
  end
  for j=1:N
    ampo += -2*h,"Sx",j;
  end
  H = MPO(ampo,sites)


  #
  # Make sure to do lots of sweeps
  # when finding excited states
  #
  sweeps = Sweeps(30)
  maxdim!(sweeps,10,10,10,20,20,40,80,100,200,200)
  cutoff!(sweeps,1E-8)
  noise!(sweeps,1E-6)

  #
  # Compute the ground state psi0
  #
  psi0_init = randomMPS(sites,2)
  energy0,psi0 = dmrg(H,psi0_init,sweeps)

  println()

  #
  # Compute the first excited state psi1
  # (Use ground state psi0 as initial state
  #  and as a 'penalty state')
  #
  psi1_init = psi0
  energy1,psi1 = dmrg(H,[psi0],psi1_init,sweeps; weight)

  # Check psi1 is orthogonal to psi0
  @show inner(psi1,psi0)


  #
  # The expected gap of the transverse field Ising
  # model is given by Eg = 2*|h-1|
  #
  # (The DMRG gap will have finite-size corrections.)
  #
  println("DMRG energy gap = ",energy1-energy0);
  println("Theoretical gap = ",2*abs(h-1));

  println()

  #
  # Compute the second excited state psi2
  # (Use ground state psi0 as initial state
  #  and [psi0,psi1] as 'penalty states')
  #
  psi2_init = psi0
  energy2,psi2 = dmrg(H,[psi0,psi1],psi2_init,sweeps;weight)

  @show inner(psi2,psi0)
  @show inner(psi2,psi1)

  return
end

# Results of the code with N = 20, h = 4:
# 1 energy=-81.183755787555 maxlinkdim=10 maxerr=1.07E-07 time=36.364
# After sweep 2 energy=-81.191708707619 maxlinkdim=3 maxerr=9.71E-09 time=0.065
# After sweep 3 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.052
# After sweep 4 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.044
# After sweep 5 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.046
# After sweep 6 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.046
# After sweep 7 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.070
# After sweep 8 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.060
# After sweep 9 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.064
# After sweep 10 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.054
# After sweep 11 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.060
# After sweep 12 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.077
# After sweep 13 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.087
# After sweep 14 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.049
# After sweep 15 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.054
# After sweep 16 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.056
# After sweep 17 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.074
# After sweep 18 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.048
# After sweep 19 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.062
# After sweep 20 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.055
# After sweep 21 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.056
# After sweep 22 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.070
# After sweep 23 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.054
# After sweep 24 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.063
# After sweep 25 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.061
# After sweep 26 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.088
# After sweep 27 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.045
# After sweep 28 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.052
# After sweep 29 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.055
# After sweep 30 energy=-81.191708714622 maxlinkdim=3 maxerr=4.86E-09 time=0.054

# After sweep 1 energy=-74.839796777601 maxlinkdim=10 maxerr=1.21E-08 time=2.617
# After sweep 2 energy=-75.001801990735 maxlinkdim=10 maxerr=1.32E-08 time=0.068
# After sweep 3 energy=-75.066935934139 maxlinkdim=10 maxerr=4.86E-08 time=0.094
# After sweep 4 energy=-75.102466692902 maxlinkdim=12 maxerr=9.72E-09 time=0.148
# After sweep 5 energy=-75.123456370723 maxlinkdim=11 maxerr=9.40E-09 time=0.085
# After sweep 6 energy=-75.139349154230 maxlinkdim=11 maxerr=9.46E-09 time=0.096
# After sweep 7 energy=-75.148984622850 maxlinkdim=10 maxerr=9.88E-09 time=0.066
# After sweep 8 energy=-75.154478422023 maxlinkdim=10 maxerr=9.69E-09 time=0.078
# After sweep 9 energy=-75.158076260928 maxlinkdim=9 maxerr=9.41E-09 time=0.085
# After sweep 10 energy=-75.160155077682 maxlinkdim=9 maxerr=8.81E-09 time=0.078
# After sweep 11 energy=-75.161281352147 maxlinkdim=8 maxerr=9.55E-09 time=0.074
# After sweep 12 energy=-75.162269704163 maxlinkdim=8 maxerr=9.77E-09 time=0.084
# After sweep 13 energy=-75.162783935068 maxlinkdim=9 maxerr=8.68E-09 time=0.068
# After sweep 14 energy=-75.162874108780 maxlinkdim=6 maxerr=9.46E-09 time=0.069
# After sweep 15 energy=-75.162895279096 maxlinkdim=6 maxerr=9.78E-09 time=0.103
# After sweep 16 energy=-75.162902838479 maxlinkdim=6 maxerr=9.99E-09 time=0.073
# After sweep 17 energy=-75.162905845965 maxlinkdim=6 maxerr=9.13E-09 time=0.070
# After sweep 18 energy=-75.162906746083 maxlinkdim=6 maxerr=9.97E-09 time=0.087
# After sweep 19 energy=-75.162907284319 maxlinkdim=6 maxerr=9.93E-09 time=0.059
# After sweep 20 energy=-75.162907379192 maxlinkdim=6 maxerr=9.88E-09 time=0.077
# After sweep 21 energy=-75.162907422283 maxlinkdim=6 maxerr=9.86E-09 time=0.077
# After sweep 22 energy=-75.162907443677 maxlinkdim=6 maxerr=9.85E-09 time=0.084
# After sweep 23 energy=-75.162907454769 maxlinkdim=6 maxerr=9.84E-09 time=0.067
# After sweep 24 energy=-75.162907460678 maxlinkdim=6 maxerr=9.84E-09 time=0.070
# After sweep 25 energy=-75.162907463968 maxlinkdim=6 maxerr=9.84E-09 time=0.086
# After sweep 26 energy=-75.162907465900 maxlinkdim=6 maxerr=9.84E-09 time=0.063
# After sweep 27 energy=-75.162907467091 maxlinkdim=6 maxerr=9.84E-09 time=0.066
# After sweep 28 energy=-75.162907467850 maxlinkdim=6 maxerr=9.84E-09 time=0.091
# After sweep 29 energy=-75.162907468347 maxlinkdim=6 maxerr=9.84E-09 time=0.080
# After sweep 30 energy=-75.162907468676 maxlinkdim=6 maxerr=9.84E-09 time=0.063
# inner(psi1, psi0) = -9.409140133698202e-15
# DMRG energy gap = 6.028801245945857
# Theoretical gap = 6.0

# After sweep 1 energy=-74.657644078277 maxlinkdim=10 maxerr=1.40E-07 time=0.055
# After sweep 2 energy=-74.864147188274 maxlinkdim=10 maxerr=1.02E-07 time=0.077
# After sweep 3 energy=-74.961948681070 maxlinkdim=10 maxerr=5.43E-08 time=0.098
# After sweep 4 energy=-75.014370584207 maxlinkdim=13 maxerr=9.54E-09 time=0.074
# After sweep 5 energy=-75.038743396976 maxlinkdim=14 maxerr=9.77E-09 time=0.108
# After sweep 6 energy=-75.053296836022 maxlinkdim=11 maxerr=9.79E-09 time=0.078
# After sweep 7 energy=-75.063906682228 maxlinkdim=11 maxerr=9.64E-09 time=0.091
# After sweep 8 energy=-75.069036929847 maxlinkdim=11 maxerr=9.54E-09 time=0.106
# After sweep 9 energy=-75.072126415937 maxlinkdim=9 maxerr=9.75E-09 time=0.077
# After sweep 10 energy=-75.074900871670 maxlinkdim=8 maxerr=9.95E-09 time=0.096
# After sweep 11 energy=-75.076069222808 maxlinkdim=8 maxerr=8.00E-09 time=0.072
# After sweep 12 energy=-75.076811272757 maxlinkdim=8 maxerr=9.25E-09 time=0.086
# After sweep 13 energy=-75.077398968295 maxlinkdim=7 maxerr=9.71E-09 time=0.087
# After sweep 14 energy=-75.077710202358 maxlinkdim=6 maxerr=9.71E-09 time=0.072
# After sweep 15 energy=-75.077794020445 maxlinkdim=6 maxerr=9.51E-09 time=0.099
# After sweep 16 energy=-75.077827385744 maxlinkdim=6 maxerr=9.39E-09 time=0.075
# After sweep 17 energy=-75.077841345364 maxlinkdim=6 maxerr=9.99E-09 time=0.138
# After sweep 18 energy=-75.077847162752 maxlinkdim=6 maxerr=9.93E-09 time=0.132
# After sweep 19 energy=-75.077849572618 maxlinkdim=6 maxerr=9.83E-09 time=0.114
# After sweep 20 energy=-75.077850527459 maxlinkdim=6 maxerr=9.95E-09 time=0.106
# After sweep 21 energy=-75.077850917301 maxlinkdim=6 maxerr=9.77E-09 time=0.121
# After sweep 22 energy=-75.077851084742 maxlinkdim=6 maxerr=9.75E-09 time=0.094
# After sweep 23 energy=-75.077851159318 maxlinkdim=6 maxerr=9.74E-09 time=0.097
# After sweep 24 energy=-75.077851193347 maxlinkdim=6 maxerr=9.74E-09 time=0.117
# After sweep 25 energy=-75.077851209061 maxlinkdim=6 maxerr=9.73E-09 time=0.099
# After sweep 26 energy=-75.077851216443 maxlinkdim=6 maxerr=9.73E-09 time=0.124
# After sweep 27 energy=-75.077851220182 maxlinkdim=6 maxerr=9.73E-09 time=0.098
# After sweep 28 energy=-75.077851222265 maxlinkdim=6 maxerr=9.73E-09 time=0.084
# After sweep 29 energy=-75.077851223513 maxlinkdim=6 maxerr=9.72E-09 time=0.137
# After sweep 30 energy=-75.077851224298 maxlinkdim=6 maxerr=9.72E-09 time=0.081
# inner(psi2, psi0) = -1.2281842209915794e-14
# inner(psi2, psi1) = -1.3558464248635514e-6
