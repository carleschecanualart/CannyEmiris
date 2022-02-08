# CannyEmiris
Greedy implementation of the Canny-Emiris formula as it was presented in "A Greedy Approach to the Canny-Emiris Formula" with the introduction of type functions. 

## Zonotopes 

The matrix A has to represent the bounds of the supports 

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{A}'_i = \big\{ \sum_{j = 1}^n \lambda_j v_j \in \mathbb{Z}^n\, | \quad \lambda_j \in \mathbb{Z}, \quad 0  \leq \lambda_j \leq a_{ij}\big\}">
while the columns of the matrix H represent the line segments defining the zonotope. The program will compute also the exponent <img src="https://render.githubusercontent.com/render/math?math=det(H)"> that appears in <img src="https://render.githubusercontent.com/render/math?math=\Res_{\mathcal{A}'} = \Res_{\mathcal{A}}^{|det(H)|}"> where 

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{A}_i = \big\{(b_j)_{j = 1,\dots,n} \in \mathbb{Z}^n \quad | \quad 0 \leq b_j \leq a_{ij} \big\} \quad i = 0,\dots,n"> 

CannyEmiris.Zonotopes(A::Matrix,H::Matrix) writes the rows of the Canny-Emiris matrix and returns two symbolic matrices which are the Canny-Emiris matrix and its principal minor.

``` console

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
julia> D
2×2 Matrix{SymPy.Sym}:
 (u_{2, [0, 0]})  (u_{2, [1, 1]})
 (u_{1, [0, 0]})  (u_{1, [1, 1]})
````

The program will return an error if the matrices do not shape to the system or if H does not correspond to an n-zonotope. It will also return an error if the entries of A are not ordered.
 	
## Multihomogeneous systems 
	
Let now N represent the vector of n_1,...n_s in a multihomogeneous system in <img src="https://render.githubusercontent.com/render/math?math=\mathbb{P}^{n_1} \times \dots \times \mathbb{P}^{n_s}"> and let D be a matrix whose columns are the multidegrees of the polynomials of the system. 
	
CannyEmiris.Zonotopes(D::Matrix,N::Vector) writes the rows of the Canny-Emiris matrix and returns two symbolic matrices which are the Canny-Emiris matrix and its principal minor.
	
```julia

julia> D
1×3 Matrix{Int64}:
 2  2  1

julia> N
1-element Vector{Int64}:
 2

julia> CE,PM = CannyEmiris.Multihomogeneous(D,N)
The rows of the Canny-Emiris matrix x^{b-a(b)}F_{i(b)} are:
[2, 1]-> x^[2, 1]*F_2
[3, 1]-> x^[3, 1]*F_2
[1, 2]-> x^[1, 2]*F_2
[2, 2]-> x^[2, 2]*F_2
[1, 3]-> x^[1, 3]*F_2
[4, 1]-> x^[2, 1]*F_1
[3, 2]-> x^[1, 2]*F_1
[2, 3]-> x^[2, 1]*F_0
[1, 4]-> x^[1, 2]*F_0

The size of the greedy Canny-Emiris matrix is: 9
The degree of the resultant is: 8

The exponent of difference between this resultant and the one given by the canonical basis of Z^n is: 1.0

julia> CE
9×9 Matrix{SymPy.Sym}:
 (u_{2, [0, 0]})  (u_{2, [1, 0]})                0  (u_{2, [0, 1]})                0                0                0                0                0
               0  (u_{2, [0, 0]})                0                0                0  (u_{2, [1, 0]})  (u_{2, [0, 1]})                0                0
               0                0  (u_{2, [0, 0]})  (u_{2, [1, 0]})  (u_{2, [0, 1]})                0                0                0                0
               0                0                0  (u_{2, [0, 0]})                0                0  (u_{2, [1, 0]})  (u_{2, [0, 1]})                0
               0                0                0                0  (u_{2, [0, 0]})                0                0  (u_{2, [1, 0]})  (u_{2, [0, 1]})
 (u_{1, [0, 0]})  (u_{1, [1, 0]})                0  (u_{1, [0, 1]})                0  (u_{1, [2, 0]})  (u_{1, [1, 1]})  (u_{1, [0, 2]})                0
               0                0  (u_{1, [0, 0]})  (u_{1, [1, 0]})  (u_{1, [0, 1]})                0  (u_{1, [2, 0]})  (u_{1, [1, 1]})  (u_{1, [0, 2]})
 (u_{0, [0, 0]})  (u_{0, [1, 0]})                0  (u_{0, [0, 1]})                0  (u_{0, [2, 0]})  (u_{0, [1, 1]})  (u_{0, [0, 2]})                0
               0                0  (u_{0, [0, 0]})  (u_{0, [1, 0]})  (u_{0, [0, 1]})                0  (u_{0, [2, 0]})  (u_{0, [1, 1]})  (u_{0, [0, 2]})

julia> PM
1×1 Matrix{SymPy.Sym}:
 (u_{2, [0, 0]})
 
 ````
 
 Here after
	
	
	
 
