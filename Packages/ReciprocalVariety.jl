# Load packages
using HomotopyContinuation, LinearAlgebra

# Set HC to don't use compilation since the computations are all small
set_default_compile(:none)



# Built index matrix for catalecticant space Cat(r,s)

sym_to_vec(L,m) = [L[i, j] for i = 1:m for j = i:m]


function vec_to_sym(l::AbstractVector{T},m::Int64) where T
    k = 0
    L = Matrix{T}(undef, m, m)
    for i in 1:m, j in i:m
        k += 1
        L[i, j] = L[j, i] = l[k]
    end
    L
end

function CatalecticantIndices(r,s,m)
    tf = [sum(I) == r for I in Iterators.product(ntuple(i->0:r, s)...)]
    indx = collect.(collect(Iterators.product(ntuple(i->0:r, s)...))[tf])
    M = hcat([indx for _ in 1:m]...)
    return M+permutedims(M)
end

function CatalecticantSpace(M::AbstractMatrix,A::AbstractMatrix,m::Integer)
    eq = Vector{typeof(A[1,1]+A[1,1])}(undef,0)
    count = Vector{typeof(M[1,1])}(undef,0) 
    ## Parsing the upper triangular matrix
    for i in 1:m
        for j in i:m
            ## Comparing with the remaining entries
            for i2 in (i+1):m
                for j2 in i2:m
                    if M[i,j] == M[i2,j2]
                       #forcing it's not compared to some entry already used to compare
                       if ([i,j] in count)==false
                         #push!(eq, A[i,j]-A[i2,j2]);
                         push!(eq, A[j,i]-A[j2,i2]);
                         push!(count,[i,j]);
                       end
                    end
                end
            end
        end
    end
return eq
end

function adjugate(M::AbstractMatrix{T}) where T
    nr, nc = size(M)
    out = similar(M)
    rows = BitArray(ones(Int16,nr))
    cols = BitArray(ones(Int16,nc))
    for r in 1:nr
        for c in r:nc
            rows[r] = 0
            cols[c] = 0
            out[c, r] = (-1)^(c+r)*det(M[rows,cols])
            rows[r] = 1
            cols[c] = 1
        end
    end
    return out
end

