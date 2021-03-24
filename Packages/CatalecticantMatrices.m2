newPackage(
        "CatalecticantMatrices",
        Version => "1.0", 
        Date => "6 December 2020",
        Authors => {
            {Name => "Laura Brustenga i Moncusi", 
                Email => "brust@math.ku.dk", 
             Name => "Elisa Cazzador", 
                Email => "elisacaz@math.uio.no", 
             Name => "Roser Homs", 
                Email => "roser.homs@tum.de", 		
		}
            },
        Headline => "tools for catalecticant matrices"
        )

export {    
    "adjugate",
    "genericCatalecticantMatrix"
    }

----------------------------------------------------------------------------------------
-- The present Macaulay2 file is an auxiliar package which is part of the             --
-- complementary material to the paper "Inverting catalecticants of ternary quartics" --                       
-- by Laura Brustenga i Moncusí, Elisa Cazzador and Roser Homs                        --
----------------------------------------------------------------------------------------


----------------------------------------------------------------------------------------
-- ADJUGATE MATRIX                                                                    --
----------------------------------------------------------------------------------------   
-- Compute adjugate of a symmetric matrix                                             --
-- Input:                                                                             --
--   * M, a matrix                                                                    --
-- Output:                                                                            --
--   * adjugate matrix of M                                                           --
----------------------------------------------------------------------------------------   

  adjugate = M -> ( 
    m := numcols M; 
    if m != numrows M then error "Matrix is not square";
    if entries M != entries transpose M then error "Matrix is not symmetric";
    adjugateM := for i to m-1 list ( 
	  for j to m-1 list (-1)^(i+j)*det(submatrix'(M, {j}, {i})) 
	  ); 
    matrix adjugateM
  ) 
 

----------------------------------------------------------------------------------------
-- GENERIC CATALECTICANT MATRIX                                                       --
----------------------------------------------------------------------------------------
-- Generic catalecticant matrix of (n+1)-ary forms of degree 2k with                  --
-- coefficients in the ring X (i.e. the entries of the matrix are the variables in X) --
-- Input:                                                                             --
--    * k, an integer                                                                      --
--    * n, an integer                                                                      --
--    * X, a ring  with indexed variables                                                 --
-- Output:                                                                            --
--    * generic catalecticant matrix with entries in the ring X                       --
----------------------------------------------------------------------------------------

  genericCatalecticantMatrix = (k,n,X) -> (
    --projective dimension of space of catalecticants
    N:=#(gens X)-1;
    --local variables for auxiliary ring R
    a:=symbol a;
    x:=symbol x;
    --ring with variables a (coefficients of the form) and x (variables of the form)
    R:=QQ[a_0..a_N,x_0..x_n];
    -- list of variables a
    aList:=drop(gens R,-(n+1));
    -- list of variables x
    xList:=drop(gens R,N+1);
    -- basis of monomials of degree d=2k in variables x
    dBasis:=flatten entries basis(2*k,R,Variables=>xList);
    -- expression of a (n+1)-ary form of degree 2k 
    -- with variables x and coefficients a
    form:=sum(for i to #aList-1 list aList_i*dBasis_i);
    -- basis of monomials of degree k in variables x
    kBasis:=flatten entries basis(k,R,Variables=>xList);
    -- contraction (apolar action) of the form 
    --with respect to operators given by monomials of degree k
    CM:=contract(matrix{kBasis},form);
    -- retrieve coefficients a of the contraction: catalecticant matrix
    (M,C):=coefficients(CM,Variables=>xList,Monomials=>kBasis);
    -- define ring map to return catalectican matrix to X
    f:=map(X,R,join(gens X,apply(xList,i->0)));
    -- retrieve catalecticant matrix in X, remove unnecessary information 
    matrix entries f(C)
  )
