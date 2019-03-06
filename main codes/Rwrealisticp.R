k = 3
d = 100
B = 500
b = 3.6
n = 1000
n0 = 25
n1 = n-n0
M = rep(0,B)
y = matrix(rnorm(d*n),n,d)
for (u in 1:B){
  y <- y[sample(nrow(y)),]
  y.dist = as.matrix(dist(y))
  diag(y.dist) = max(y.dist)+100
  
  An = matrix(0,n,k)
  for (i in 1:n){
    An[i,] = (sort(y.dist[i,], index.return=T)$ix)[1:k]
  }
  
  temp = table(An)
  id = as.numeric(row.names(temp))
  deg = rep(0,n)
  deg[id] = temp
  cn = sum((deg-k)^2)/n/k
  R2 = (cn + k)*n*k
  count = 0
  for (i in 1:n){
    ids = An[i,]
    count = count + length(which(An[ids,]==i))
  }
  vn = count/n/k
  
  Z = rep(0,n)
  
  for (t in n0:n1){
    EX = -(k*n*(t - 1)*(t - n + 1))/(n^2 - 3*n + 2)
    VX = (t*(t - 1)*(n^2 - 2*n*t - n + t^2 + t)*(R2 + 2*count - R2*n - 3*count*n + 2*k*n + count*n^2 - 3*k*n^2 + 3*k^2*n + k*n^3 - k^2*n^2))/(n*(n - 3)*(n^2 - 3*n + 2)^2)
    a1 = (n-t-1)/(n-2)
    a2 = (t-1)/(n-2)
    Z[t] = (a1*length(which(An[1:t,]<=t))+a2*length(which(An[(t+1):n,]>t))-EX)/sqrt(VX)
  }
  M[u] = max(Z)}

prob = length(which(M>=b))/B

