----------------------------------------------------------------------------------------
-- The present Macaulay2 file contains complementary material                         --
-- to the paper "Inverting catalecticants of ternary quartics"                        --
-- by Laura Brustenga i Moncusí, Elisa Cazzador and Roser Homs                        --
--                                                                                    --
-- Complementary material for Proposition 3.5, Remark 4.4 and Proposition 4.5         --
----------------------------------------------------------------------------------------

restart
loadPackage "CatalecticantMatrices" 
-- The package CatalecticantMatrices.m2 is required for the following functions:
-- genericCatalecticantMatrix, adjugate


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

-- Build ring of the parametrization, including the variable t
Z = QQ[x_0..x_N,y_0..y_M,t] 

cat = sub(cat,Z)

-- Representative for the open orbit of rank-4 points,
-- lying on the 4-secant to v4(1:0:0), v4(0:1:0), v4(0:0:1) and v4(1:1:1)
rank4 = sub(cat, join( 
	{x_0=>2, x_10=>2, x_N=>2}, 
	for i from 1 to 9 list x_i=>1,
	for i from 11 to N-1 list x_i=>1  
    ))

-- Representative for the open orbit of rank-3 points
-- lying on the 4-secant to v4(1:0:0), v4(0:1:0) and v4(0:0:1)
rank3 = sub(cat, join(
	{x_0=>1, x_10=>1, x_N=>1}, 
	for i from 1 to 9 list x_i=>0,
	for i from 11 to N-1 list x_i=>0  
    ))

-- Representative for the open orbit of rank-2 points on a proper secant
-- lying on the 2-secant to v4(1:0:0) and v4(0:0:1)
rank2Secant = sub(cat, join(
    {x_0=>1, x_N=>1}, 
    for i from 1 to N-1 list x_i=>0 
    ))

-- Representative for the closed orbit of rank-2 points on a tangent line 
-- lying on a tangent line to v4(1:0:0)
rank2Tangent = sub(cat, join(
    {x_0=>1, x_1=>1},
    for i from 2 to N list x_i=>0
    ))

-- Representative for the open orbit of rank-1 points 
-- v4(1:0:0)
rank1 = sub(cat, join(
    {x_0=>1},
    for i from 1 to N list x_i=>0
    ))

orbits = {rank4, rank3, rank2Secant, rank2Tangent, rank1};

-- Build the list of parametrizations for each type of matrix A listed above.
-- The parametrization is deduced from LEMMA 3.4
parametrizations = {}; 
for A in orbits do(
    r := rank A;
    line := A + t*cat;
    imageLine := adjugate line; -- phi(A + t*Cat)
    divImageLine := matrix(
	for i to m-1 list (for j to m-1 list imageLine_(i,j) // 
	    t^(m-r-1))
	); -- divide all the entries by the maximum power of t dividing them
    limit := sub(divImageLine, t=>0); -- set t=0
    entriesLimit := flatten(
	for i to m-1 list(
	    for j from i to m-1 list limit_(i,j)
	    )
	);
    parametrizePhiA := ideal(
	for i to M list y_i - entriesLimit_i
	); -- parametrization of phi(A)
    parametrizations = parametrizations | {parametrizePhiA};
    )

-- Build the list containing the defining equations for the image phi(A) for r>1
equations = {}; 
for i to 3 do(    
        phiA := sub(
	    eliminate(toList(x_0..x_N), parametrizations_i), 
	    Y); -- takes about 10 minutes -- elimination to get equations
        equations = equations | {phiA};
	)

-- Print projective dimensions, degree and defining ideals of phi(A) when r > 1
for E in equations do print (dim E - 1, degree E, E)

-- Compare the two classes of cubic 8-folds with the cubic Pfaffians of the two
-- skew-symmetric matrices S1 and S2
secantCubic = ideal (equations_2)_11
tangentCubic = ideal (equations_3)_11

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

-- When r=1, there are 6 easy linear equations in the parametrization of phi(A) 
use Z
paramRank1 = parametrizations_4
PP14 = ideal(for i to 5 list paramRank1_i) -- 6 linear equations 

-- When r=1, the dimension of phi(A) is 11. We check this by intersecting the
-- parametrization of phi(A) with the pull-back of 11 random hyperplanes of PS^m
use Z
H11 = ideal flatten entries( -- equations of 11 random hyperplanes
    random(QQ^11, QQ^21)*transpose(matrix{toList(y_0..y_M)}) 
    );
intersec = paramRank1 + H11;
elim = sub(
    eliminate(toList(x_0..x_N), intersec), 
    Y); -- takes approx. 30 minutes
dim elim-1, degree elim


----------------------------------------------------------------------------------------
-- REMARK 4.4: A point in the orthogonal of PCat(2,3) has generically full rank.      --
-- The orthogonal space is empty in rank 1 and 2, a Veronese surface in rank 3 and    --
-- a cubic hypersurface, secant to the Veronese, in rank 4 and 5.                     --
----------------------------------------------------------------------------------------

-- Build "ring of PCat(2,3) x PS^m", where the y_i' are seen as coefficients
XY = Y[x_0..x_N] 
cat = sub(cat, XY);
sym = sub(sym, XY);

-- Build orthogonal space to PCat(2,3) and its rank loci
traceProduct = trace (cat*sym)
orthoCat = trim ideal(
    for i to N list(
	sub(coefficient(x_i, traceProduct), Y)
	)
    ) 
dim orthoCat-1

sym = sub(sym, Y);
ortho1 = orthoCat + minors(2, sym);
ortho2 = orthoCat + minors(3, sym);
ortho3 = orthoCat + minors(4, sym);
ortho4 = orthoCat + minors(5, sym);
ortho5 = orthoCat + ideal(det sym);

-- A general point in the orthogonal haas full rank
rank sub(sym, Y/orthoCat)

-- The orthogonal space is empty in rank 1 and 2
dim ortho1-1, dim ortho2-1

-- Build Veronese surface and its secant and compare with loci of rank 3,4,5
use Y
mat = matrix{
    {y_18, y_13, y_12},
    {y_13, y_11, y_7},
    {y_12, y_7, y_6}
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
sym = sub(sym, Y);
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
