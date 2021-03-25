# Load file for needed functions
include("../Packages/ReciprocalVariety.jl")


# Basic setup. For ternary quartics: r=2, s=3
    r=2 
    s=3
    m=binomial(r+s-1,s-1)
    N=binomial(2r+s-1,2r)
    M=binomial(m + 1, 2)
    @var x[0:2r, 0:2r, 0:2r],y[1:M]

# Define ideal equations
    cat=CatalecticantIndices(r,s,m)
    adj=adjugate(Expression.(vec_to_sym(y,m)));
    eq=CatalecticantSpace(cat,adj,m);

# Set up system
    F = System(eq; variables=y);

# Compute point on the variety
    x₀ = randn(ComplexF64,N)
    @var x[0:2r, 0:2r, 0:2r]
    X = map(ijk -> x[(ijk .+ 1)...], cat)
    y₀ = sym_to_vec(inv(X(variables(X) => x₀)),m)
    
# Put affine subspace through this point
    V₀ = let
       A₀ = randn(ComplexF64, N, M)
       b₀ = A₀ * y₀
       LinearSubspace(A₀, b₀)
    end

# Perform monodromy
    mres = monodromy_solve(
       F,
       y₀,
       V₀,
       parameter_sampler=_ -> LinearSubspace(
           randn(ComplexF64, N, M),
           randn(ComplexF64, N),
       ),
   )



