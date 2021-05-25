using Yao
c= yao"""let nqubits=9, version="0.6.0"
    begin # encode circuit
        1=>C, 4=>X
        1=>C, 7=>X
        1=>H, 4=>H, 7=>H
        1=>C, 2=>X
        1=>C, 3=>X
        4=>C, 5=>X
        4=>C, 6=>X
        7=>C, 8=>X
        7=>C, 9=>X
    end

    # the error
    1=>X, 2=>Z, 3=>Z, 4=>X, 5=>Z, 6=>Z, 7=>X, 8=>Z, 9=>Z

    being # decode circuit
        1=>C, 2=>X
        1=>C, 3=>X
        2=>C, 3=>C, 1=>X
        4=>C, 5=>X
        4=>C, 6=>X
        5=>C, 6=>C, 4=>X
        7=>C, 8=>X
        7=>C, 9=>X
        8=>C, 9=>C, 7=>X
        1=>H, 4=>H, 7=>H
        1=>C, 4=>X
        1=>C, 7=>X
        4=>C, 7=>C, 1=>X
    end
end"""
