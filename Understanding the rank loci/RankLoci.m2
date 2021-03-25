----------------------------------------------------------------------------------------
-- The present Macaulay2 file contains complementary material                         --
-- to the paper "Inverting catalecticants of ternary quartics"                        --
-- by Laura Brustenga i Moncusí, Elisa Cazzador and Roser Homs                        --
--                                                                                    --
-- Complementary material for Proposition 3.5, Remark 4.4 and Proposition 4.5         --
----------------------------------------------------------------------------------------

restart
load "../Packages/CatalecticantMatrices.m2" 

-- The package CatalecticantMatrices.m2 is required for the following functions:
-- genericCatalecticantMatrix, adjugate, parametrizedImage, toCat


----------------------------------------------------------------------------------------
-- SETUP COMPUTATIONS: CATALECTICANT AND SYMMETRIC MATRICES                           --
----------------------------------------------------------------------------------------

-- Choose space of catalecticants of (n+1)-ary forms of degree 2k
-- n=2 and k=2 corresponds to catalecticants associated with ternary quartics Cat(2,3)
n = 2  
k = 2 

-- Build generic catalecticant matrix in PCat(k,n+1) = P^N
N = binomial(n+2*k,2*k)-1                    -- projective dimension of catalecticant space
X = QQ[x_0..x_N]                           -- ring of catalecticant matrices
cat=genericCatalecticantMatrix(k,n,X)      -- catalecticant matrix

-- Build generic symmetric matrix in PS^m = P^M
m = binomial(k+n,k)                        -- size of catalecticant matrix
M = binomial(m+1,2)-1                      -- projective dimension of symmetric space
Y = QQ[y_0..y_M]                           -- ring of symmetric matrices
sym = genericSymmetricMatrix(Y,m)          -- symmetric matrix


----------------------------------------------------------------------------------------
-- PROPOSITION 3.5: The image of a rank-r point is a P^2, a P^5, a Pfaffian cubic     --
-- 8-fold in P^9 and an 11-fold in P^14 for r=4,3,2,1 respectively. We prove this     --
-- using a parametrization deduced from LEMMA 3.4                                     --
----------------------------------------------------------------------------------------

-- STEP 1 
-- Parametrize the image of points in the rank-r open orbits for r=1..4 and in the 
-- closed orbit of rank-2 points                                               

-- Points in P^2 to build the proper secants
P={1,0,0}, Q={0,1,0}, R={0,0,1}, S={1,1,1}

-- Representatives for the orbits
rank4 = toCat(P,4) + toCat(Q,4) + toCat(R,4) + toCat(S,4)
rank3 = toCat(P,4) + toCat(Q,4) + toCat(R,4)
rank2Secant = toCat(P,4) + toCat(R,4)
rank2Tangent = sub(cat, join({x_0=>1, x_1=>1}, for i from 2 to N list x_i=>0 ))
rank1 = toCat(P,4)

orbits = {rank4, rank3, rank2Secant, rank2Tangent, rank1};

-- List of parametrizations for every orbit representative
Z = X**Y;
parametrizations = for A in orbits list(
    sub(parametrizedImage(A, cat, sym), Z)
    );


-- STEP 2
-- Use elimination to find equations for the point images when r>1
-- and compute dimension and degree of the resulting varieties

use Z
equations =  for i to 3 list(
    sub(eliminate(toList(x_0..x_N), parametrizations_i), Y) -- takes about 10 minutes
    );

for E in equations do print (dim E - 1, degree E, E)


-- STEP 3
-- Check that the two classes of cubic 8-folds are the cubic Pfaffian of a 6x6 skew-
-- symmetric matrix

-- Two classes of cubics
secantCubic = ideal (equations_2)_11
tangentCubic = ideal (equations_3)_11

-- Build skew-symmetric matrices
use Y
S1 = matrix{
    {0, 0, y_6, y_7, y_8, y_9},
    {0, 0, y_7, y_11, y_12, y_13},
    {-y_6, -y_7, 0, y_12-y_9, y_15, y_16},
    {-y_7, -y_11, y_9-y_12, 0, y_16, y_18},
    {-y_8, -y_12, -y_15, -y_16, 0, 0},
    {-y_9, -y_13, -y_16, -y_18, 0, 0}
    }
S2 = matrix{
    {0, 0, y_19, y_20, y_14, y_17},
    {0, 0, y_13, y_14, y_11, y_12},
    {-y_19, -y_13, 0, y_17-y_18, y_12, y_15},
    {-y_20, -y_14, y_18-y_17, 0, y_13, y_16},
    {-y_14, -y_11, -y_12, -y_13, 0, 0},
    {-y_17, -y_12, -y_15, -y_16, 0, 0}
    }

pfaffians(6, S1) == secantCubic
pfaffians(6, S2) == tangentCubic


-- STEP 4
-- Case r=1. Check that the image of a rank-1 point  is an 11-fold embedded in a P^14
 
-- We easily find 6 linear equations in the parametrization
use Z
paramRank1 = parametrizations_4
PP14 = ideal(for i to 5 list paramRank1_i) 

-- Intersect the parametrization with the pull-back of 11 random hyperplanes of PS^m
use Z
pullback11 = ideal flatten entries( -- equations of 11 random hyperplanes
    random(QQ^11, QQ^21)*transpose(matrix{toList(y_0..y_M)}) 
    );
intersec = paramRank1 + pullback11;

-- After elimination, get a zero-dimensional variety
elim = sub(
    eliminate(toList(x_0..x_N), intersec), 
    Y); -- takes approx. 30 minutes
dim elim-1, degree elim


----------------------------------------------------------------------------------------
-- REMARK 4.4: A point in the orthogonal of PCat(2,3) has generically full rank.      --
-- The orthogonal space is empty in rank 1 and 2, a Veronese surface in rank 3 and    --
-- a cubic hypersurface, secant to the Veronese, in rank 4 and 5.                     --
----------------------------------------------------------------------------------------

-- Build the orthogonal space to PCat(2,3) and its rank loci
orthoCat = orthogonal(cat, sym)
dim orthoCat-1

rk1ortho = orthoCat + minors(2, sym);
rk2ortho = orthoCat + minors(3, sym);
rk3ortho = orthoCat + minors(4, sym);
rk4ortho = orthoCat + minors(5, sym);
rk5ortho = orthoCat + ideal(det sym);

-- A general point in the orthogonal haas full rank
rank sub(sym, Y/orthoCat)

-- The orthogonal space is empty in rank 1 and 2
dim ortho1-1, dim ortho2-1

-- Build a suitable Veronese surface and its secant and compare with loci of rank 3,4,5
use Y
mat = matrix{
    {y_6, y_7, y_12},
    {y_7, y_11, y_13},
    {y_12, y_13, y_18}
    }
veronese = minors(2, mat)
secant = ideal det(mat)

radical ortho3 == veronese + orthoCat
radical ortho4 == secant + orthoCat
radical ortho5 == secant + orthoCat


---------------------------------------------------------------------------------------
-- PROPOSITION 4.5: The orthogonal of PCat(2,3) does not contain any full-rank point --
-- and intersects phi(C_r) in the emptyset for r=3,4,5 and a Veronese surface r=1,2  --
---------------------------------------------------------------------------------------

-- Ideal of the pull-back phi^*(PCat(2,3)) in PS^m
use Y;
pullbackCat = trim ideal(
    det submatrix'(sym, {0}, {3}) + det submatrix'(sym, {1}, {1}),
    det submatrix'(sym, {0}, {4}) + det submatrix'(sym, {1}, {2}), 
    det submatrix'(sym, {0}, {5}) + det submatrix'(sym, {2}, {2}),  
    det submatrix'(sym, {1}, {4}) - det submatrix'(sym, {2}, {3}),  
    det submatrix'(sym, {1}, {5}) - det submatrix'(sym, {2}, {4}),  
    det submatrix'(sym, {3}, {5}) - det submatrix'(sym, {4}, {4})  
    );

-- Ideal of the intersection of phi^*(PCat(2,3)) with the orthogonal space
intersectPreimageCat = trim(pullbackCat + orthoCat);

-- Symmetric matrices in such intersection
intersectionMatrix = sub(sym, Y/intersectPreimageCat) 

-- Matrices in the intersection are degenerate
det intersectionMatrix == 0 

-- Intersection with phi(A), rk(A) = 2,3
PP5 = equations_1;
secant8fold = equations_2;
tangent8fold = equations_3;

dim (PP5+orthoCat)-1
dim (secant8fold+orthoCat)-1
dim (tangent8fold+orthoCat)-1, degree (tangent8fold+orthoCat)
