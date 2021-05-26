using DifferentialEquations

#Removed some code that didn't work
#Make another attempt if I find time
"f(t,y) = 0.5y
u0 = 1.5
timespan = (0.0, 1.0) #Solve from time = 0 to time = 1
prob = ODEProblem(f,u0,timespan)

sol = solve(prob) # Solves the ODE

A = [1. 0 0 -5
    4 -2 4 -3
    -4 0 0 1
    5 -2 2 3]
u0 = rand(4,2)
f(t,u) = A*u0
prob = ODEProblem(f,u0,timespan)
sol = solve(prob)"


# The famous Lorenz system in a more friendly way
f = @ode_def Lorenz begin
    dx = σ*(y-x)
    dy = x*(ρ - z) - y
    dz = x*y - β*z
end σ ρ β
# needed to change from σ=>10. ρ=>28. β=(8/3)
