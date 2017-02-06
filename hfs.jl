"""Computes harmonic function solution"""
function iterative_hfs(A::SparseMatrixCSC,Y_masked::Vector{Int64},l::Int64,T::Int64;tol = 0.50)

  num_samples = size(A,1);
  num_classes = sum(unique(Y_masked) .!= 0);
  B = sparse(exp(-A) .* (A .> 0));
  sumB = sum(B,1);

  # compute the target y for the linear system
  l_idx = find(Y_masked);
  u_idx = find(Y_masked.==0);

  y = zeros(Int64,length(l_idx),num_classes);
  for k=1:num_classes
    y[:,k] = Y_masked[l_idx].==k;
  end

  # compute the hfs solution, using iterated averaging
  f = zeros(num_samples,num_classes);
  f[l_idx,:] = y;

  for t=1:T
    for k=1:num_classes
      f[u_idx,k] = (sum(spdiagm(f[:,k],0,num_samples,num_samples)*B[:,u_idx],1)./sumB[u_idx]')';
    end
  end

  #= compute the labels assignment from the hfs solution
  label: (n x 1) class assignments [0,1,2,...,num_classes] =#
  labels = zeros(Int64,num_samples);
  labels[l_idx] = Y_masked[l_idx];
  m,idx = findmax(f[u_idx,:],2);
  labels[u_idx] = ind2sub(size(f[u_idx,:]), vec(idx))[2] .* (m .> tol);


  return labels

end

# Mask labels from a vector
# l number of unmasked (revealed) labels to include in the output
# r number of communities
function mask_labels(Y::Vector{Int64}, l::Int64, r::Int64)

  num_samples = size(Y,1);
  nb = div(num_samples,r);

  Y_masked = copy(Y);

  for i = 1:r
    idx = randperm(nb);
    Y_masked[(i-1)*nb+idx[l+1:end]] = 0;
  end

  return Y_masked
end
