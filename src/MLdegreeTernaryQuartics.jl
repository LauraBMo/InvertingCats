
include("CatalecticantMatrices.jl")

### Create the model
Σ = Catalecticant(2, 3)

### Compute the ML-degree of Σ
ml_degree_witness(Σ; dual=true)
# MLDegreeWitness:
#  • ML degree → 36
#  • model dimension → 15
#  • dual → true
