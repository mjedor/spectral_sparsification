"""calculate effective resistances in a given graph"""
function effectiveResistances(B::SparseMatrixCSC{Int64,Int64},a::Vector{Float64};epsilon::Float64=0.1)

  # Creating data for effective resitances
  m,n = size(B);
  W = spdiagm(a); # Weight matrix
  L = B'*W*B; # Laplacian matrix
  #scale_ = ceil(Int64,log2(n)/epsilon);
  scale_ = ceil(Int64,24*log(n)/epsilon^2)
  delta = epsilon/3 * sqrt(2(1-epsilon)*minimum(a)/((1+epsilon)*size(B,2)^3*maximum(a)));

  # Finding the effective resitances
  Q = (rand(scale_,m)) .> 0.5;
  Q = Q - !Q;
  Q = Q./sqrt(scale_);
  Y = sparse(Q*sqrt(W)*B); # Create the system
  W = 0; # Clear W
  Q = 0; # Clear Q

  # Solving the systems
  #= SDDSolvers = - augTreeSolver
                  - KMPSDDSolver
                  - samplingSDDSolver
                  - AMGSolver
  =#
  excess = spzeros(n); excess[1] = excess[n] = 0.1;
  f = KMPSDDSolver(L+spdiagm(excess),tol=delta);
  Z = spzeros(scale_,size(L)[1])
  L = 0; # Clear L

  for j=1:scale_
    Z[j,:] = f(full(Y[j,:]));
  end
  Y = 0; # Clear Y

  # Creating list of edges
  elist = zeros(Int64,size(B,1),2);
  for i=1:size(B,1)
    elist[i,:] = [find(B[i,:].==-1) find(B[i,:].==1)]
  end

  return sum(((Z[:,elist[:,1]]-Z[:,elist[:,2]]).^2),1)'
end

function effectiveResistances(H::Sparsifier;epsilon=0.01)
  return effectiveResistances(H.B,H.a,epsilon=epsilon)
end
