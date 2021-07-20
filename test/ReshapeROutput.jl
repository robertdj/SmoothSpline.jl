function vector_to_banded_matrix(x)
    nrows = div(length(x), 4)

    X = zeros(nrows, nrows)

    X[LinearAlgebra.diagind(X)] = x[1:nrows]

    X[LinearAlgebra.diagind(X, 1)] = x[nrows + 1:2*nrows - 1]
    X[LinearAlgebra.diagind(X, -1)] = x[nrows + 1:2*nrows - 1]

    X[LinearAlgebra.diagind(X, 2)] = x[2*nrows + 1:3*nrows - 2]
    X[LinearAlgebra.diagind(X, -2)] = x[2*nrows + 1:3*nrows - 2]

    X[LinearAlgebra.diagind(X, 3)] = x[3*nrows + 1:4*nrows - 3]
    X[LinearAlgebra.diagind(X, -3)] = x[3*nrows + 1:4*nrows - 3]

    return X
end


function vector_to_upper_triangular(x)
    nrows = div(length(x), 4)

    X = zeros(nrows, nrows)

    # factorization in R with dpbfa.f. Storage format is peculiar; works backwards
    X[1, 1] = x[4]
    X[1:2, 2] = x[7:8]
    X[1:3, 3] = x[10:12]
    for col in 4:nrows
        idx1 = col - 4 .+ (1:4)
        idx2 = (col - 1) * 4 .+ (1:4)
        X[idx1, col] = x[idx2]
    end

    return X
end

