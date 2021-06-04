using ITensors

#=
This code is written to provide the energy levels, and state psi if desired,
of an arbitrary excited state. Note: If it's a highly excited state, it'll
have to go through all states below it to accumulate the "penalty states" so
the code is inefficient for highly excited states.
=#

function excited_state(N, h, n)
  # N denotes number of qubits, h is proportional to the "weight",
  # and n is the number of iterations to go through.

  # It returns the energy level of the nth excited state.
  sites = siteinds("S=1/2",N)

  let
    weight = 10*h
    ampo = AutoMPO()
    for j=1:N-1
      ampo += -4,"Sz",j,"Sz",j+1
    end
    for j=1:N
      ampo += -2*h,"Sx",j;
    end
    H = MPO(ampo,sites)

    sweeps = Sweeps(20)
    maxdim!(sweeps,10,10,10,20,20,40,80,100,140,90, 20)
    cutoff!(sweeps,1E-8)
    noise!(sweeps,1E-6)

    psi0_init = randomMPS(sites,2)
    println("Ground state energy is:")
    energy0,psi0 = dmrg(H,psi0_init,sweeps)
    println()

    penalty_states = [psi0]
    energies = [energy0]
    for i = 1:n
      psii_init = penalty_states[1]
      println("Energy level = ", i, " is:")
      energyi, psii = dmrg(H, penalty_states, psii_init, sweeps; weight)
      println()
      push!(penalty_states, psii)
      push!(energies, energyi)
    end

    if string(n) == string(11) || string(n) == string(12) || string(n) == string(13)
      println(n,"th excited energy level is ",energies[end])
    elseif string(n)[end] == Char('1')
      println(n,"st excited energy level is ",energies[end])
    elseif string(n)[end] == Char('2')
      println(n,"nd excited energy level is ",energies[end])
    elseif string(n)[end] == Char('3')
      println(n,"rd excited energy level is ",energies[end])
    else
      println(n,"th excited energy level is ",energies[end])
    end
    println()
    # N.B. The code doesn't correctly call the excited states 11th, 12th, 13th
    # as they don't follow the usual pattern; extra conditions would be needed
    # to correctly name these states, but it's not essential to how the code
    # runs so this hasn't been done for simplicity.

  end
end

excited_state(25, 5, 12)

# This code produces the same result with excited_state(20, 4, 2) as
# the ITensor excited states code produces, but is now generalised
# to any number of states.
