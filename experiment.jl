using Laplacians
using Distributions
using DataStructures

include("project.jl")
using project

# Seed
srand(123)

# Parameters for graph construction
r = 3; # number of communities
nb = 40; # number of vertices per communities
n = repmat([nb],r);
P = rand(Beta(1,3),r,r);
for i = 1:r
    P[i,i] = rand(Beta(5,1),1,1)[1];
end

# Parameters for sparsification
k = 8;
epsilon = 1.;
alpha = 1 + epsilon;
delta = 0.1;
N=sum(n);
qBar = ceil(Int64,12*alpha*log(2*alpha*N/delta)/epsilon^2);

# Parameters for hfs
T = 7;
l = 4;
tol = 0.20;

# Construct label vector
Y = zeros(Int64,N);
nb = div(N,r);
for i = 1:r
  Y[(i-1)*nb+1:i*nb] = i;
end

# Mask some labels
Y_masked = mask_labels(Y,l,r);

# Construct SBM graph
A = stochastic_block_model(P,n);
A = betaWeight(A,n);
isConnected(A)

# Sparsify
@time H = sparsify(A,k,qBar,epsilon=epsilon)

# Compute hfs
avg = mean(Y .== iterative_hfs(A,Y_masked,l,T,tol = tol))

# Compute hfs with sparsifier
Lh = lapFromEdge(H.B,H.a);
Ah,Dh = decomposeLaplacian(Lh);
avgS = mean(Y .== iterative_hfs(Ah,Y_masked,l,T,tol = tol))
