"""Create an adjacency matrix and a diagonal vector from a Laplacian"""
function decomposeLaplacian{Tv,Ti}(la::SparseMatrixCSC{Tv,Ti})
  a = -(2 * la - triu(la) - tril(la))
  d = diag(la)
  return a,d
end

"""Set the weight of every edge to random uniform [0,1]"""
function uniformWeightSym!(A::SparseMatrixCSC)
  (ai,aj) = findnz(triu(A,1))
  for i in 1:length(ai)
      A[ai[i],aj[i]] = rand(1)[1]
      A[aj[i],ai[i]] = A[ai[i],aj[i]]
  end
end

"""The signed edge-vertex adjacency matrix and weight vector"""
function edgeVertexMatWithWeight(A::SparseMatrixCSC)
    (ai,aj,a) = findnz(triu(A,1))
    m = length(ai)
    n = size(A)[1]
    return convert(SparseMatrixCSC{Int64,Int64},sparse(collect(1:m),ai,1.0,m,n) - sparse(collect(1:m),aj,1.0,m,n)),a
end

function lapFromEdge(B::SparseMatrixCSC,a::Vector{Float64})
  return B' * spdiagm(a) * B;
end
