<<<<<<< HEAD
using Yao

N = 6
Xstring = chain(N, prod([put(N, i=>X) for i=1:N]))
ZZpairs = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])

iscommute(Xstring, ZZpairs)

U(Jt::Float64) = time_evolve(0.5 * Xstring + 0.5 * ZZpairs, Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function Mz_evolve(nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U(i * deltaJt)
        else
            ψ |> U(i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

using Plots
t_vec, Mz_vec = Mz_evolve(10000, 0.00001)
plot(t_vec[9000:10000], Mz_vec[9000:10000])
=======
using Yao

N = 6
Xstring = chain(N, prod([put(N, i=>X) for i=1:N]))
ZZpairs = sum([chain(N, put(N, i=>Z)*put(N, i+1=>Z)) for i=1:N-1])

U(Jt::Float64) = time_evolve(0.5 * Xstring + 0.5 * ZZpairs, Jt, tol=1e-5, check_hermicity=true)
Mz(N::Int) = sum([put(N, i => Z) for i = 1:N]) / N

function Mz_evolve(nsteps, deltaJt)
    protected = true
    t_vec = Vector{Float64}();
    Mz_vec = Vector{Float64}();
    ψ = zero_state(N)
    for i = 0:nsteps
        append!(t_vec, i * deltaJt)
        if protected
            ψ |> U(i * deltaJt)
        else
            ψ |> U(i * deltaJt)
        end
        append!(Mz_vec, expect(Mz(N), ψ))
    end
    return t_vec, Mz_vec
end

using Plots
t_vec, Mz_vec = Mz_evolve(10000, 0.00001)
plot(t_vec[9000:10000], Mz_vec[9000:10000])
>>>>>>> 53ebef02eaf0b1d0ce3196697ed00f740ff38d99
