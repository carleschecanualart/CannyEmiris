# CannyEmiris
Greedy implementation of the Canny-Emiris formula as it was presented in "A Greedy Approach to the Canny-Emiris Formula" with the introduction of type functions. The underlying system of equations is written as:

<img src="https://render.githubusercontent.com/render/math?math=F_i = \sum_{a \in \mathcal{A}_i}u_{i,a}\chi^a \quad i = 0,\dots,n">

for some finite sets of supports <img src="https://render.githubusercontent.com/render/math?math=\mathcal{A}_i \subset M"> in a lattice of rank n. This implementation treats the cases in which the Newton polytopes 
<img src="https://render.githubusercontent.com/render/math?math=\Delta_i = conv(\mathcal{A}_i)"> are zonotopes (sums of line segments) or products of simplices (which correspond to multihomogeneous systems of equations).

## Zonotopes 

Let <img src="https://render.githubusercontent.com/render/math?math=v_1,\dots,v_n \in M"> be independent vectors generating an n-zonotope. The supports of our zonotopes are: 

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{A}_i^' = \big\{ \sum_{j = 1}^n \lambda_j v_j \in \mathbb{Z}^n\, | \quad \lambda_j \in \mathbb{Z}, \quad 0  \leq \lambda_j \leq a_{ij}\big\}">

for some <img src="https://render.githubusercontent.com/render/math?math=a_{i,j} \quad i = 0,\dots,n \quad j = 1,\dots,s">. Let V be the <img src="https://render.githubusercontent.com/render/math?math=n \times n"> matrix whose columns are <img src="https://render.githubusercontent.com/render/math?math=v_1,\dots,v_n"> and let A be the matrix whose entries are the <img src="https://render.githubusercontent.com/render/math?math=a_{i,j}">.

The function CannyEmiris.Zonotopes receives these two matrices and returns two symbolic matrices <img src="https://render.githubusercontent.com/render/math?math=\mathcal{H}_{\mathcal{G}},\mathcal{E}_{\mathcal{G}}"> which correspond to the rational formula for the sparse resultant:

<img src="https://render.githubusercontent.com/render/math?math=Res_{\mathcal{A}} = \big(\frac{det(\mathcal{H}_{\mathcal{G}})}{det(\mathcal{E}_{\mathcal{G}})}\big)^{det(V)}">

which follows from Theorem 1.1 and Corollary 3.1 on the text. 

Moreover, the program prints the polynomials <img src="https://render.githubusercontent.com/render/math?math=\chi^{b-a(b)}F_{i(b)}"> corresponding to all the lattice points in the greedy subset <img src="https://render.githubusercontent.com/render/math?math=b \in \mathcal{G}">, the size of the matrix and the degree of the resultant (which corresponds to the lattice points in mixed cells).

Let's see an example of the use of this function which corresponds to the system of Example 1.1 in the text.

``` julia

julia> include("CannyEmiris.jl")
Main.CannyEmiris

julia> A = [[1,1] [1,1] [1,1]]
2×3 Matrix{Int64}:
 1  1  1
 1  1  1

julia> H = [[1,0] [0,1]]
2×2 Matrix{Int64}:
 1  0
 0  1

julia> CE, PM = CannyEmiris.Zonotopes(A,H)
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

The sparse resultant is the ratio of the determinants of the returned matrices to the power 1.0


julia> CE
8×8 Matrix{SymPy.Sym}:
 (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})                0                0                0
 (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})                0                0                0
               0                0  (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})                0
               0                0                0  (u_{2, [0, 0]})  (u_{2, [0, 1]})                0  (u_{2, [1, 0]})  (u_{2, [1, 1]})
 (u_{0, [0, 0]})  (u_{0, [0, 1]})                0  (u_{0, [1, 0]})  (u_{0, [1, 1]})                0                0                0
               0                0  (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})                0
               0                0  (u_{0, [0, 0]})  (u_{0, [0, 1]})                0  (u_{0, [1, 0]})  (u_{0, [1, 1]})                0
               0                0                0  (u_{1, [0, 0]})  (u_{1, [0, 1]})                0  (u_{1, [1, 0]})  (u_{1, [1, 1]})

julia> PM
2×2 Matrix{SymPy.Sym}:
 (u_{2, [0, 0]})  (u_{2, [1, 1]})
 (u_{1, [0, 0]})  (u_{1, [1, 1]})
````

The program will return an error if i) the matrix A or H do not have the desired dimensions, 
``` julia
julia> A = [[1,1] [1,1]]
2×2 Matrix{Int64}:
 1  1
 1  1

julia> H
2×2 Matrix{Int64}:
 1  0
 0  1

julia> CE, PM = CannyEmiris.Zonotopes(A,H)
The matrix of the a_{i,j} does not have the correct dimensions
(Dict{Any, Any}(), Dict{Any, Any}())
````
ii) if H does not correspond to an n-zonotope

``` julia
julia> H = [[1,0] [0,0]]
2×2 Matrix{Int64}:
 1  0
 0  0

julia> CE, PM = CannyEmiris.Zonotopes(A,H)
The vectors do not correspond to an n-zonotope.
(Dict{Any, Any}(), Dict{Any, Any}())
```

or iii) entries of A are not ordered in the sense  <img src="https://render.githubusercontent.com/render/math?math=0 < a_{0j} \leq a_{1j} \leq \dots \leq a_{n-1,j} \quad j = 1,\dots,n.">:

```julia
julia> A = [[2,1] [1,1] [1,1]]
2×3 Matrix{Int64}:
 2  1  1
 1  1  1

julia> H = [[1,0] [0,1]]
2×2 Matrix{Int64}:
 1  0
 0  1

julia> CE, PM = CannyEmiris.Zonotopes(A,H)
The matrix of the a_{i,j} does not satisfy a_{i-1,j} <= a_{i,j}
(Dict{Any, Any}(), Dict{Any, Any}())
```
 	
## Multihomogeneous systems 
	
Let now N represent the vector of n_1,...n_s in a multihomogeneous system in <img src="https://render.githubusercontent.com/render/math?math=\mathbb{P}^{n_1} \times \dots \times \mathbb{P}^{n_s}"> and let D be a matrix whose columns are the multidegrees of the polynomials of the system so that:

<img src="https://render.githubusercontent.com/render/math?math=\mathcal{A}_i = \big\{(b_{jl})_{l = 1,\dots,s}^{j = 1,\dots,n_l} \in \oplus_{l=1}^s\mathbb{Z}^{n_l} | b_{jl} \geq 0 \quad \sum_{j=0}^{n_l}b_{jl} \leq d_{i,l} \} \quad i = 0,\dots,n">
	
The function CannyEmiris.Multihomogeneous(D::Matrix,N::Vector) receives these two matrices D and N and returns and returns two symbolic matrices which are the Canny-Emiris matrix and its principal minor, as before. Let's show its use in the system corresponding to Example 4.1.
	
```julia
julia> N = [2]
1-element Vector{Int64}:
 2

julia> D = [2 2 1]
1×3 Matrix{Int64}:
 2  2  1

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

(SymPy.Sym[(u_{2, [0, 0]}) (u_{2, [1, 0]}) … 0 0; 0 (u_{2, [0, 0]}) … 0 0; … ; (u_{0, [0, 0]}) (u_{0, [1, 0]}) … (u_{0, [0, 2]}) 0; 0 0 … (u_{0, [1, 1]}) (u_{0, [0, 2]})], SymPy.Sym[(u_{2, 
[0, 0]});;])

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
 
 Again, this function will return nothing if the dimension conditions and the order. 
 
## Other functions to call

· MultihomogeneousEmbedding
	
The function CannyEmiris.MultihomogeneousEmbedding(MULTI_A::Matrix{Int64}, MULTI_N::Vector{Int64})) builds the embedding of the multihomogeneous system into a zonotope system. This system is given by the supports:

<img src="https://render.githubusercontent.com/render/math?math=\overline{\mathcal{A}_i} = \big\{(b_{jl})_{l = 1,\dots,s}^{j=1,\dots,n_l} \in \oplus_{j=1}^s\mathbb{Z}^{n_j} \: | \: 0 \leq \sum_{j = J}^{n_l}b_{jl} \leq d_{i,j}\}\quad l = 1,\dots,s \quad J = 1,\dots,n_l">

```julia

julia> N = [2]
1-element Vector{Int64}:
 2

julia> D = [2 2 1]
1×3 Matrix{Int64}:
 2  2  1

julia> CannyEmiris.MultihomogeneousEmbedding(D,N)
([2 2 1; 2 2 1], [1 -1; 0 1])
 
 ````
 
 · GenerateTypeFunctions
 
The function CannyEmiris.GenerateTypeFunctions(n::Int) constructs the iterator that produces all the type functions <img src="https://render.githubusercontent.com/render/math?math=\varphi:\{1,\dots,n\} \xrightarrow{} \{0,\dots,n\}"> satisfying the condition on the greedy subset <img src="https://render.githubusercontent.com/render/math?math=\mathcal{G}">:

<img src="https://render.githubusercontent.com/render/math?math=\sum_{i = 0}^{I-1} t_{b,i} \leq I \quad \forall I = 1,\dots,n"> where <img src="https://render.githubusercontent.com/render/math?math=t_{b,i} = |\varphi^{-1}(i)|">

```julia

julia> p = CannyEmiris.GenerateTypeFunctions(2)
Main.CannyEmiris.TypeFunctions{Vector{Int64}}([0, 1, 2], 2)

julia> for x in p println(x) end
[0, 1]
[0, 2]
[1, 0]
[1, 1]
[1, 2]
[2, 0]
[2, 1]
[2, 2]
 
 ````
