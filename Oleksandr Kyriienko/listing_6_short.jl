using Yao
using YaoExtensions
using KrylovKit: eigsolve;

h = heisenberg(10)
w, v = eigsolve(mat(h), 10, :SR, ishermitian=true)

println("Eigenenergies are $w")
