# Type Sparsifier
type Sparsifier
  B::SparseMatrixCSC{Int64,Int64}
  a::Vector{Float64}
  p::Vector{Float64}
  q::Vector{Int64}
end

Sparsifier(B::SparseMatrixCSC{Int64,Int64},a::Vector{Float64},q::Int64) = Sparsifier(B,a,ones(length(a)),q.*ones(length(a)))
