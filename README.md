# CannyEmiris
Greedy implementation of the Canny-Emiris formula as it was presented in "A Greedy Approach to the Canny-Emiris Formula" with the introduction of type functions.

CannyEmiris.Zonotopes(A::Matrix,H::Matrix) writes the rows of the Canny-Emiris matrix and 

julia> include("CannyEmiris.jl")

julia> A = [[1,1] [1,1] [1,1]]
2×3 Matrix{Int64}:
 1  1  1
 1  1  1

julia> H = [[1,0] [0,1]]
2×2 Matrix{Int64}:
 1  0
 0  1
 
julia> C, D = CannyEmiris.Zonotopes(A,H)
The rows of the Canny-Emiris matrix x^{b-a(b)}F_{i(b)} are:
[0, 1]-> x^[0, 1]*F_2
[0, 2]-> x^[0, 1]*F_1
[1, 0]-> x^[1, 0]*F_2
[1, 1]-> x^[1, 1]*F_2
[1, 2]-> x^[0, 1]*F_0
[2, 0]-> x^[1, 0]*F_1
[2, 1]-> x^[1, 0]*F_0
[2, 2]-> x^[1, 1]*F_1

The size of the greedy Canny-Emiris matrix is: 8
The degree of the resultant is: 6

The exponent of difference between this resultant and the one given by the canonical basis of Z^n is: 1.0

julia> C
8×8 Matrix{SymPy.Sym}:
 (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})                0                0                0
 (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})                0                0                0
               0                0  (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})                0
               0                0                0  (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})
 (u_{0, [0, 0]})  (u_{0, [0, 1]})                0  (u_{0, [1, 0]})  (u_{0, [1, 1]})                0                0                0
               0                0  (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})                0
               0                0  (u_{0, [0, 0]})  (u_{0, [0, 1]})                0  (u_{0, [1, 0]})  (u_{0, [1, 1]})                0
               0                0                0  (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})
2×2 Matrix{SymPy.Sym}:
 (u_{2, [0, 0]})  (u_{2, [1, 1]})
 (u_{1, [0, 0]})  (u_{1, [1, 1]})             
 
