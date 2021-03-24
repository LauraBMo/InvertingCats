----------------------------------------------------------------------------------------
-- The present Macaulay2 file contains complementary material                         --
-- to the paper "Inverting catalecticants of ternary quartics"                        --
-- by Laura Brustenga i Moncusí, Elisa Cazzador and Roser Homs                        --
--                                                                                    --
-- Complementary material for the proofs of Lemma 3.4 and Proposition 4.4             --
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

-- Build generic catalecticant matrix in in P^N, namely PCat(k,n+1) 
N = binomial(n+2*k,2*k)-1                    -- projective dimension of catalecticant space
X = QQ[x_0..x_N]                           -- ring of catalecticant matrices
cat=genericCatalecticantMatrix(k,n,X)      -- catalecticant matrix

-- Build generic symmetric matrix of the same size in P^M
m = binomial(k+n,k)                        -- size of catalecticant matrix
M = binomial(m+1,2)-1                      -- projective dimension of symmetric space
Y = QQ[y_0..y_M]                           -- ring of symmetric matrices
sym = genericSymmetricMatrix(Y,m)          -- symmetric matrix


----------------------------------------------------------------------------------------
-- LEMMA 3.4: The dimension of phi(C2) is 12. We prove it by showing that for every   --
-- matrix A of rank 2 we have dim( pi_2(pi_1^(-1)(A) \cap Gamma ) = 8                 --
----------------------------------------------------------------------------------------

-- Build ring of the parametrization, including the variable t
Z = QQ[x_0..x_N,y_0..y_M,t] 

cat = sub(cat,Z)

-- Representative for the orbit of rank-2 point corresponding to a proper secant
type1 = sub(cat, 
    {x_0=>1, x_N=>1} | 
    for i from 1 to N-1 list x_i=>0 
    )

-- Representative for the orbit of rank-2 point corresponding to a tangent line
type2 = sub(cat,
    {x_0=>1, x_1=>1} |
    for i from 2 to N list x_i=>0
    )

orbits = {type1, type2};

-- Build the list containing pi_2(pi_1^(-1)(A) \cap Gamma) 
-- for each rype of rank-2 matrix A 
adjOrbits = {};
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
    parametrizeAdjA := ideal(
	for i to M list y_i - entriesLimit_i
	); -- parametrization of pi_2(pi_1^(-1)(A) \cap Gamma)
    adjA := sub(
	eliminate(toList(x_0..x_N), parametrizeAdjA), 
	Y); -- takes about 10 minutes -- get expicit equations via elimination
    adjOrbits = adjOrbits | {adjA};
    )

-- Print projective dimensions, degree and defining ideals of the two cubic 8-folds
for B in adjOrbits do print (dim B - 1, degree B, B)


----------------------------------------------------------------------------------------
-- PROPOSITION 4.4: PCat(2,3)^(-1) intersects the ortogonal only on points of phi(C2) --
-- that are associated with tangent lines to C1 (and possibly along phi(C1))          --
----------------------------------------------------------------------------------------

-- Build "ring of P^N x P^M", where the y_i' are seen as coefficients
XY = Y[x_0..x_N] 
cat = sub(cat, XY);
sym = sub(sym, XY);

-- Build orthogonal space to the P^N of catalecticant matrices
traceProduct = trace (cat*sym)
orthoCat = ideal(
    for i to N list(
	sub(coefficient(x_i, traceProduct), Y)
	)
    ) 
dim orthoCat - 1


---------------------------------------------------------------------------------------
-- Step 1: Check that for any full rank matrix A in P^N, its image phi(A) is not in  --
-- the orthogonal space i.e. the preimage phi^(-1)(P^14) intersects the orthogonal   --
-- in degenerate matrices                                                            --
---------------------------------------------------------------------------------------

-- Ideal of the preimage phi^(-1)(P^14) in P^20
preimageCat = trim ideal(
    det submatrix'(sym, {0}, {3}) - det submatrix'(sym, {1}, {1}),
    det submatrix'(sym, {0}, {4}) - det submatrix'(sym, {1}, {2}), 
    det submatrix'(sym, {0}, {5}) - det submatrix'(sym, {2}, {2}),  
    det submatrix'(sym, {1}, {4}) - det submatrix'(sym, {2}, {3}),  
    det submatrix'(sym, {1}, {5}) - det submatrix'(sym, {2}, {4}),  
    det submatrix'(sym, {3}, {5}) - det submatrix'(sym, {4}, {4})  
    );
preimageCat = sub(preimageCat, Y);

-- Ideal of the intersection of phi^(-1)(P^14) with the orthogonal space
intersectPreimageCat = trim(preimageCat + orthoCat);

-- Symmetric matrices in such intersection
intersectionMatrix = sub(sym, Y/intersectPreimageCat) 

-- Matrices in the intersection are degenerate
det intersectionMatrix == 0 


--------------------------------------------------------------------------------------
-- Step 2: Check how the two cubic 8-folds of Lemma 3.4  intersect the orthogonal   -- 
-- space                                                                            --
--------------------------------------------------------------------------------------

-- Print projective dimension and degree the two cubic 8-folds intersected with the 
-- orthogonal space
for B in adjOrbits do(
    intersection := B + orthoCat;
    print (dim intersection - 1, degree intersection)
    )

