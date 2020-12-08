# InvertingCats
Repository containing the computations of the paper *Inverting catalecticant matrices of ternary quartics* by Laura Brustenga i Moncusí, Elisa Cazzador and Roser Homs.

The code we present here involves two essential steps to understand the relation between the degree of the reciprocal variety of Cat(2,3) and the ML-degree of the model represented by Cat(2,3).

**Computing the ML-degree**
(*Section 4*)

Julia code in [MLDegreeTernaryQuartics.jl](https://github.com/LauraBMo/InvertingCats/blob/main/src/MLdegreeTernaryQuartics.jl), using the package [LinearCovarianceModels.jl](https://github.com/saschatimme/LinearCovarianceModels.jl), provides an efficient procedure to compute the ML-degree of Cat(2,3). This value is 36. In fact, this can be used to compute the ML-degree of any space Cat(k,n+1) of catalecticant matrices of (n+1)-ary forms of degree 2k.

The Macaulay2 code in [MLdegree.m2](https://github.com/LauraBMo/InvertingCats/blob/main/src/MLdegree.m2) presents a symbolic aproach to the same calculation. However, it only provides an actual answer for spaces of catalecticant matrices associated to binary forms of small degree.

**Understanding the locus of rank 2 matrices**
(*Lemma 3.4 and Proposition 4.4*)

Macaulay2 computations in [Rank2Locus.m2](https://github.com/LauraBMo/InvertingCats/blob/main/src/Rank2Locus.m2) allow us to prove that the dimension of the rank 2 locus of the inversion map is exactly 12 (Lemma 3.4) and that it intersects the orthogonal of Cat(2,3) (Proposition 4.4). Therefore, the degree of the reciprocal variety of Cat(2,3) has 36 as a strict lower bound.
