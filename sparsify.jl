function sparsify(A::SparseMatrixCSC,k::Int64,qBar::Int64;epsilon=0.1)
  # Partion graph
  B_l,a_l = partitionGraph(A,k);

  # Construct initial set of sparsifiers
  S = Queue(Sparsifier);
  for l=1:k
    enqueue!(S,Sparsifier(B_l[l],a_l[l],qBar));
  end

  # Clear B_l, a_l
  B_l = 0;
  a_l = 0;

  # Merge tree
  for h=1:k-1
    H1 = dequeue!(S);
    H2 = dequeue!(S);
    enqueue!(S,mergeResparsify(H1,H2,qBar,epsilon=epsilon))
    gc();
  end

  return dequeue!(S)

end
