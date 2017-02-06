"""Merge-Resparsify function"""
function mergeResparsify(H1::Sparsifier,H2::Sparsifier,qBar::Int64;epsilon=0.1)

  # Compute effectives resistances
  r = effectiveResistances([H1.B;H2.B],[H1.a;H2.a],epsilon=epsilon);

  # Compute new weights
  a = ([H1.a;H2.a] .* [H1.q;H2.q]) ./ (qBar * [H1.p;H2.p]);

  # Compute new probability
  p = vec(min([H1.a;H2.a].*r,[H1.p;H2.p]));

  # Merge H1 and H2 (H = H1 + H2)
  H = Sparsifier([H1.B;H2.B],a,[H1.p;H2.p],[H1.q;H2.q]);

  # Resparsify
  e=1;
  while e <= size(H.B,1)
    # Sample q_e from a binomial
    Bin = Binomial(H.q[e],p[e]/H.p[e]);
    q_e = rand(Bin,1)[1];

    if(q_e==0)
      # Drop edge from H
      H.B = H.B[[1:e-1;e+1:end;],:];
      splice!(H.a,e);
      splice!(H.p,e);
      splice!(H.q,e);
      splice!(p,e);
    else
      # Update edge
      H.p[e] = p[e];
      H.q[e] = q_e;
      e+=1;
    end
  end

  # Clear memory
  H1 = 0;
  H2 = 0;
  #gc();

  return H
end
