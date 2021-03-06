# InvertingCats
Repository containing the computations of the paper *Inverting catalecticant matrices of ternary quartics* by Laura Brustenga i Moncusí, Elisa Cazzador and Roser Homs.

The code we present here involves two essential steps to understand the relation between the degree of the reciprocal variety of Cat(2,3) and the ML-degree of the model represented by Cat(2,3).

**Computing the ML-degree**
(*Section 2*)

Julia code in [MLDegreeTernaryQuartics.jl](https://github.com/LauraBMo/InvertingCats/blob/main/Computing%20the%20ML-degree/MLdegreeTernaryQuartics.jl), using the package [LinearCovarianceModels.jl](https://github.com/saschatimme/LinearCovarianceModels.jl), provides an efficient procedure to compute the ML-degree of Cat(2,3). This value is 36. In fact, this can be used to compute the ML-degree of any space Cat(k,n+1) of catalecticant matrices of (n+1)-ary forms of degree 2k.

The Macaulay2 code in [MLdegree.m2](https://github.com/LauraBMo/InvertingCats/blob/main/Computing%20the%20ML-degree/MLdegree.m2) presents a symbolic aproach to the same calculation. However, it only provides an actual answer for spaces of catalecticant matrices associated to binary forms of small degree (see Example 1.5 of the paper).

**Computing the degree**
(*Section 2*)

Homotopy continuation techniques allow us to compute the degree of the reciprocal varieties of catalecticant matrices associated to binary forms of small degree and ternary quartics in [DegreeTernaryQuartics.jl](https://github.com/LauraBMo/InvertingCats/blob/main/Computing%20the%20degree/DegreeTernaryQuartics.jl). We use the Julia package [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) to that purpose. 


**Understanding the rank loci**
(*Sections 3 and 4*)

Macaulay2 computations in [RankLoci.m2](https://github.com/LauraBMo/InvertingCats/blob/main/Understanding%20the%20rank%20loci/RankLoci.m2) allow us to describe the image of a rank-r catalecticant matrix for any r=1,...,5 (Proposition 3.5) and its intersection with the orthogonal of Cat(2,3) (Proposition 4.5). In particular, since not all these intersections are empty, the degree of the reciprocal variety must be strictly greater than the ML-degree of Cat(2,3).
