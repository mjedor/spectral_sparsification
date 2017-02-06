"""Returns a Graph generated according to the Stochastic Block Model (SBM)"""
# n : number of vertices in each block
# P : matrix of probabilities inter-communities
function stochastic_block_model(P::Matrix{Float64}, n::Vector{Int}; seed::Int = -1)
  # init dsfmt generator without altering GLOBAL_RNG
  seed = -1;
  seed > 0 && Base.dSFMT.dsfmt_gv_init_by_array(MersenneTwister(seed).seed+1)
  rng = seed > 0 ? MersenneTwister(seed) : MersenneTwister()

  N = sum(n);
  A = spzeros(N,N);

  K = length(n);
  nedg = zeros(Int,K, K)

  cum = [sum(n[1:a]) for a=0:K]

  for a=1:K
    ra = cum[a]+1:cum[a+1]
    for b=a:K
      m = a==b ? div(n[a]*(n[a]-1),2) : n[a]*n[b]
      p = a==b ? n[a]*P[a,b] / (2m) : n[a]*P[a,b]/m
      nedg = rand(Binomial(m, P[a,b]))
      rb = cum[b]+1:cum[b+1]
      i=0
      while i < nedg
        source = rand(rng, ra)
        dest = rand(rng, rb)
        if source != dest
          A[source,dest] = 1;
          A[dest,source] = A[source,dest];
          i += 1;
        end
      end
    end
  end
  return A
end

function betaWeight(A::SparseMatrixCSC,n::Vector{Int})

  N = sum(n)
  K = length(n)
  cum = [sum(n[1:a]) for a=0:K]

  A1 = rand(Beta(5,1),N,N) .* A

  for i = 1:K
      A1[cum[i]+1:cum[i+1],cum[i]+1:cum[i+1]] = rand(Beta(1,3),n[i],n[i]) .* A[cum[i]+1:cum[i+1],cum[i]+1:cum[i+1]]
  end
  return A1
end
