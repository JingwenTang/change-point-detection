k = 5
d = 100
B = 1000
b = 2.7
n = 1000
n0 = 25
n1 = n-n0
M = rep(0,B)
y = matrix(rnorm(d*n),n,d)
for (u in 1:B){
  y <- y[sample(nrow(y)),]
  ##y = rbind(y[2:n,],y[1,])
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
count = 0
for (i in 1:n){
  ids = An[i,]
  count = count + length(which(An[ids,]==i))
}
vn = count/n/k

Z = rep(0,n)
for (t in n0:n1){
  EX = 2*k*t*(n-t)/(n-1)
  h = 4*(t-1)*(n-t-1)/((n-2)*(n-3))
  VX = EX*(h*(1+vn-2*k/(n-1))+(1-h)*cn)/2
  Z[t] = -(length(which(An[1:t,]>t))+length(which(An[(t+1):n,]<=t))-EX)/sqrt(VX)
}
M[u] = max(Z)}

prob = length(which(M>=b))/B
 

