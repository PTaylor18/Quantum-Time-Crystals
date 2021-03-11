using Yao
using Random
using Plots
using FFTW

Xstr(N::Int) = chain(N, prod([put(N, i=>X) for i=1:N]))
RXstr(N::Int) = chain(N, prod([put(N, i=>Rx(0.3)) for i=1:N]))
#== or ==#
XstrR(N::Int) = chain(N, repeat(X, 1:N))
Hzz(N::Int) = sum([-1 * ((put(N, i+1 => Z) * put(N, i => Z)) + put(N, i => (0.3)*Z) + put(N, i => (0.3)*X)) for i = 1:N-1])
Uzz(N::Int, Jt::Float64) = time_evolve(Hzz(N), Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function Mz_evolve(N::Int, nsteps::Int64, Jt)
    protected = true # Can change this to run the alternative XstrR and Rxstr
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i)
        if protected
            println("Hello")
            ψ |> XstrR(N) |> Uzz(N, Jt) |> RXstr(N)
        else
            ψ |> XstrR(N) |> RXstr(N)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

t_vec, Mz_vec = Mz_evolve(5, 100, π/4 #=(π)/2=#)
plot(t_vec, Mz_vec, linetype=:steppre, legend=false)
xlabel!("Time T")
ylabel!("Magnetisation Mz")

Mz_vec_fft = fft(Mz_vec)
plot(abs.(Mz_vec_fft))
