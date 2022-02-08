module CannyEmiris

export GenerateTypeFunctions, ParseMultihomogeneousEmbedding, BuildLatticeToVarZonotopes, BuildLatticeToVarMultihomogeneous, CheckMatrices
export RowsCannyEmiris, Zonotopes, CannyEmirisMultihomogeneous

using LinearAlgebra
using SymPy

struct TypeFunctions{T}
    a::T
    t::Int
end

const typefunctions(elements::T, len::Integer, rc::Integer) where {T} =
    TypeFunctions{T}(unique(elements), len)

Base.IteratorSize(::TypeFunctions) = Base.HasLength()
Base.length(p::TypeFunctions) = length(p.a)^p.t
Base.eltype(::TypeFunctions{T}) where {T} = T


function Base.iterate(p::TypeFunctions, s = [n for n = 1:p.t])


    if s[1] > lastindex(p.a)
        return nothing
    end



    cur = p.a[s]

    not_in_greedy = false

    while (not_in_greedy == false)

        i = lastindex(p.a) - 1

        s[i] += 1

        while i > 1 && s[i] > length(p.a)
            s[i] = 1
            s[i-1] += 1
            i -= 1
        end

        new_value = length(p.a) - 1  ## check new_value in G.

        while new_value > -1
            if length(findall(x -> x <= new_value, s)) >= new_value + 1
                new_value = -2
            end
            new_value -= 1
        end

        if (new_value == -1)
            not_in_greedy = true
        end

    end

    return cur, s
end


function GenerateTypeFunctions(n::Int)
    p = TypeFunctions(collect(0:n), n)
    return p
end


## Zonotopes

function BuildLatticeToVarZonotopes(A::Matrix, H::Matrix)

    n = size(A)[1]

    LatticeToVar = Dict([])

    if (size(A)[2] != n + 1 || size(H)[1] != n || size(H)[2] != n)
        return LatticeToVar = Dict([])
    end

    SymMatrix = zeros(Sym, n, n)

    for i = 1:n+1
        local base = ntuple(j -> (0:A[j, i]), n)
        v = Iterators.product(base...)
        for x in v
            coeff = symbols(
                "u_{" * string(i - 1) * "," * string(H * [j for j in x]) * "}",
            )
            var = coeff
            LatticeToVar = merge(
                LatticeToVar,
                Dict([(vcat([i - 1], H * [j for j in x]), var)]),
            )
        end
    end

    return LatticeToVar

end

function CheckMatrices(A::Matrix, H::Matrix)
    if (size(A)[1] != size(H)[1] || size(A)[1] != (size(A)[2] - 1)) return false end


    if (size(H)[1] != size(H)[2] ||det(H) == 0) return false end
    n = size(A)[1]

    for j in 1:n
        for i in 1:n-1
            if (A[j,i] > A[j,i+1]) return false
            end
        end
    end

    return true
end

function BuildLatticeToVarMultihomogeneous(A::Matrix, H::Matrix)

    n = size(A)[1]

    LatticeToVar = Dict([])

    if (size(A)[2] != n + 1 || size(H)[1] != n || size(H)[2] != n)
        return LatticeToVar
    end

    SymMatrix = zeros(Sym, n, n)

    for i = 1:n+1
        local base = ntuple(j -> (0:A[j, i]), n)
        v = Iterators.product(base...)
        for x in v
            if (all(>=(0), H * [x...]) == true)
                coeff = symbols(
                    "u_{" * string(i - 1) * "," * string(H * [x...]) * "}",
                )
                var = coeff
                LatticeToVar = merge(
                    LatticeToVar,
                    Dict([(vcat([i - 1], H * [x...]), var)]),
                )
            end
        end
    end

    return LatticeToVar
end

function BoundMatrix(A::Matrix) ## This function calculates all the partial sums \sum_{i = 0}^{I-1}a_{ij} that define the type functions

    n = size(A)[1]

    B = zeros(Int8, n, n + 2)
    for i = 1:n
        for j = 2:n+2
            B[i, j] = B[i, j-1] + A[i, j-1]
        end
    end

    return B

end

## Get the matrix H corresponding to the zonotope in which you embedded the previous system

function ParseMultihomogeneousEmbedding(MULTI_A::Matrix{Int64}, MULTI_N::Vector{Int64})

    v = ones(Int64, MULTI_N[1])
    w = -ones(Int64, MULTI_N[1] - 1)
    H = Matrix{Int64}(Bidiagonal(v, w, :U))


    for i = 2:length(MULTI_N)
        v = ones(Int64,MULTI_N[i])
        w = -ones(Int64,MULTI_N[i] - 1)
        H = cat(H, Matrix{Int64}(Bidiagonal(v, w, :U)), dims = (1, 2))
    end

    ## Get the matrix A corresponding to the degrees of the bounds of the zonotopes of H

    A = MULTI_A[1, :]'

    for i = 1:length(MULTI_N)
        for j = 1:MULTI_N[i]
            A = vcat(A, MULTI_A[i, :]')
        end
    end

    A = A[setdiff(1:end, 1), :]

    return A, H

end

## This procedure writes the rows of the Canny-Emiris matrix corresponding to G.
## At the end writes the size of the matrix with subdivision, the size of the greedy matrix and the resultant degree
## The determinant of H corresponds to the difference between our resultant and the canonical

function RowsCannyEmiris(A::Matrix, H::Matrix, ZM::Int)

    LatticeToVar = Dict([])

    number_of_rows = 0

    resultant_degree = 0

    RowToLattice = Dict([])

    RowToContent = Dict([])

    non_mixed_indices = []

    n = size(A)[1]

    p = TypeFunctions(collect(0:n), n)

    B = BoundMatrix(A)

    if (ZM == 1)
        LatticeToVar = BuildLatticeToVarZonotopes(A, H)
    end

    if (ZM == 2)
        LatticeToVar = BuildLatticeToVarMultihomogeneous(A, H)
    end

    if (length(LatticeToVar) == 0)
        return nothing
    end

    for x in p
        l = length(x)

        tb = [length(findall(x -> x == i, x)) for i = 0:l]

        non_mixed_cell = true

        ib = findlast(x -> x == 0, tb) - 1

        if (length(findall(x -> x == 1, tb)) > 1)
            non_mixed_cell = false
        end

        a = zeros(Int64, l)

        for j = 1:l
            if (x[j] > ib)
                a[j] = A[j, ib+1]
            end
        end

        local base

        if (ZM == 1)
            base = ntuple(i -> (B[i, x[i]+1]:B[i, x[i]+2]-1), n)
        end

        if (ZM == 2)
            base = ntuple(i -> (B[i, x[i]+1]+1:B[i, x[i]+2]), n)
        end

        v = Iterators.product(base...)


        for latticepoint in v
            if (ZM == 1 || (all(>=(1), H * [latticepoint...]) == true))

                number_of_rows += 1

                if (non_mixed_cell == true)
                    push!(non_mixed_indices, number_of_rows)
                end

                RowToLattice = merge(
                    RowToLattice,
                    Dict([(number_of_rows, vcat(ib, H * [latticepoint...]))]),
                )

                RowToContent =
                    merge(RowToContent, Dict([(number_of_rows, H * a)]))

                if (non_mixed_cell == false)
                    resultant_degree += 1
                end

                print(H * [latticepoint...])
                print("-> x^")
                print(H * ([latticepoint...] - a))
                print("*F_")
                println(ib)

            end

        end

    end

    println()
    print("The size of the greedy Canny-Emiris matrix is: ")
    println(number_of_rows)

    print("The degree of the resultant is: ")

    println(resultant_degree)

    println()

    print("The exponent of difference between this resultant and the one given by the canonical basis of Z^n is: ")

    println(det(H))

    return RowToLattice,
    RowToContent,
    non_mixed_indices,
    number_of_rows,
    LatticeToVar

end

function Zonotopes(A::Matrix, H::Matrix)

    if (CheckMatrices(A,H) == false)
        println("The matrices do not satisfy the bound or dimension restrictions imposed from the greedy algorithm measure.")
        return nothing
    end

    RowToLattice,
    RowToContent,
    non_mixed_indices,
    number_of_rows,
    LatticeToVar = RowsCannyEmiris(A, H, 1)

    if (RowToLattice == false) return nothing end

    n = size(A)[1]

    CannyEmirisMatrix = zeros(Sym, number_of_rows, number_of_rows)

    for i = 1:number_of_rows
        row_content = get(RowToLattice, i, 2)[1]
        row_content_support = get(RowToContent, i, 2)
        lattice_point_row = get(RowToLattice, i, 2)[2:n+1]
        for j = 1:number_of_rows
            lattice_point_column = get(RowToLattice, j, 2)[2:n+1]
            entry = Sym(
                get(
                    LatticeToVar,
                    vcat(
                        row_content,
                        lattice_point_column - lattice_point_row +
                        row_content_support,
                    ),
                    0,
                ),
            )
            CannyEmirisMatrix[i, j] = entry
        end
    end

    CannyEmirisMinor = CannyEmirisMatrix[non_mixed_indices, non_mixed_indices]

    return CannyEmirisMatrix, CannyEmirisMinor

end

function Multihomogeneous(A::Matrix, H::Vector)

    A, H = ParseMultihomogeneousEmbedding(A, H)

    if (CheckMatrices(A,H) == false)
        println("The matrices do not satisfy the bound or dimension restrictions imposed from the greedy algorithm measure.")
        return nothing
    end

    n = size(A)[1]

    RowToLattice,
    RowToContent,
    non_mixed_indices,
    number_of_rows,
    LatticeToVar = RowsCannyEmiris(A, H, 2)

    CannyEmirisMatrix = zeros(Sym, number_of_rows, number_of_rows)

    for i = 1:number_of_rows
        row_content = get(RowToLattice, i, 2)[1]
        row_content_support = get(RowToContent, i, 2)
        lattice_point_row = get(RowToLattice, i, 2)[2:n+1]
        for j = 1:number_of_rows
            lattice_point_column = get(RowToLattice, j, 2)[2:n+1]
            entry = Sym(
                get(
                    LatticeToVar,
                    vcat(
                        row_content,
                        lattice_point_column - lattice_point_row +
                        row_content_support,
                    ),
                    0,
                ),
            )
            CannyEmirisMatrix[i, j] = entry
        end
    end

    CannyEmirisMinor = CannyEmirisMatrix[non_mixed_indices, non_mixed_indices]

    return CannyEmirisMatrix, CannyEmirisMinor

end

end # module
