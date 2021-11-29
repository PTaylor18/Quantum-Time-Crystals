using IonSim
using QuantumOptics: timeevolution, stochastic
using PyPlot

C = Ca40(["S-1/2", "D-1/2"])
L = Laser()

chain = LinearChain(
        ions=[C], com_frequencies=(x=3e6,y=3e6,z=1e6),
        vibrational_modes=(;z=[1])
    )

T = Trap(configuration=chain, B=4e-4, Bhat=ẑ, δB=0, lasers=[L]);

L.k = (x̂ + ẑ)/√2
L.ϵ = (x̂ - ẑ)/√2;
Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
L.Δ = Δf

Efield_from_pi_time!(2e-6, T, 1, 1, ("S-1/2", "D-1/2"));  # Sets pi_time to 2 μs

h = hamiltonian(T, timescale=1e-6);

tspan = 0:0.01:10
mode = T.configuration.vibrational_modes.z[1]
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ionstate(T, "S-1/2") ⊗ mode[0], h);

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex, "c")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");
tight_layout()

cd(dirname(@__FILE__));
dir = pwd();
mkpath(string(dir,"\\Graphs\\"))
savefig(string(dir,"\\Graphs\\","\\Carrier_rabi_oscillations.png"))
