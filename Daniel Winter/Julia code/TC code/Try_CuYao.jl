using CuYao

r = cu(ArrayReg(bit"10101010101010"))
circ = chain(14, repeat(H, 1:14))
r |> circ
