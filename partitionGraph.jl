function partitionGraph(B::SparseMatrixCSC{Int64,Int64},a::Vector{Float64},k::Int64)

  m_l = div(size(B,1),k); # All sub-graphs have almost the same number of edges
  idx = randperm(size(B,1));
  B_l = SparseMatrixCSC{Int64,Int64}[];
  a_l = Vector{Float64}[];

  for t=1:k-1
    B_l = push!(B_l,B[idx[(t-1)*m_l+1:t*m_l],:]);
    a_l = push!(a_l,a[idx[(t-1)*m_l+1:t*m_l]]);
  end
  B_l = push!(B_l,B[idx[(k-1)*m_l+1:end],:]);
  a_l = push!(a_l,a[idx[(k-1)*m_l+1:end]]);

  return B_l,a_l

end

function partitionGraph(A::SparseMatrixCSC{Float64,Int64},k::Int64)
  B,a = edgeVertexMatWithWeight(A);
  return partitionGraph(B,a,k)
end
