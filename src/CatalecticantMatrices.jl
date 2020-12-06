

# Suit of functions to create LCModels of Catalecticant matrices.
# For example, for ternary quartics:
#
# Σ = Catalecticant(2, 3)
#
# To compute the ML-degree of Σ
#
# ml_degree_witness(Σ; dual=true)
#

using LinearCovarianceModels
import HomotopyContinuation
const HC = HomotopyContinuation

function nvectors_sumingupm(n, m)
    ntuples = Iterators.product(ntuple(_ -> 0:m, n)...)
    return filter(I -> sum(I) == m, collect(ntuples))
end

function Catalecticantindices(r, s)
    rowindeces = collect.(nvectors_sumingupm(s, r))
    m = binomial(r + s - 1, s - 1)
    M = reduce(hcat, fill(rowindeces, m))
    return M + permutedims(M)
end

function Catalecticantmatrix(r, s)
    return [HC.Variable(:θ, I...) for I in Catalecticantindices(r, s)]
end

Catalecticant = LCModel ∘ Catalecticantmatrix
