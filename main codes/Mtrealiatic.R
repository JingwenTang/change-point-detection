k = 3
d = 5
B = 500
b = 3.4
n = 1000
n0 = 200
n1 = n-n0
M = rep(0,B)
prob = rep(0,3)

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
  
  Zw = rep(0,n1-n0+1)
  Zdiff = rep(0,n1-n0+1)
  S = rep(0,n1-n0+1)
  for (t in n0:n1){
    EXw = -(k*n*(t - 1)*(t - n + 1))/(n^2 - 3*n + 2)
    VXw = (t*(t - 1)*(n^2 - 2*n*t - n + t^2 + t)*(R2 + 2*count - R2*n - 3*count*n + 2*k*n + count*n^2 - 3*k*n^2 + 3*k^2*n + k*n^3 - k^2*n^2))/(n*(n - 3)*(n^2 - 3*n + 2)^2)
    a1 = (n-t-1)/(n-2)
    a2 = (t-1)/(n-2)
    Zw[t-n0+1] = (a1*length(which(An[1:t,]<=t))+a2*length(which(An[(t+1):n,]>t))-EXw)/sqrt(VXw)
    EXdiff = -k*(n - 2*t)
    VXdiff = (t*(- n*k^2 + R2)*(n - t))/(n*(n - 1))
    Zdiff[t-n0+1] = (length(which(An[1:t,]<=t))-length(which(An[(t+1):n,]>t))-EXdiff)/sqrt(VXdiff)
    S[t-n0+1] = max(abs(Zdiff[t-n0+1]) , Zw[t-n0+1])
  }
  M[u] = max(S)}

prob[k] = length(which(M>=b))/B



