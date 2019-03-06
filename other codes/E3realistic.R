rm(list = ls())
k = 3
d = 100
b = 3.05
n = 1000
n0 = 25
n1 = n-n0
B = 100


M = matrix(0,B,n1-n0+1)
M1 = matrix(0,B,n1-n0+1)
M2 = matrix(0,B,n1-n0+1)
M3 = matrix(0,B,n1-n0+1)
M4 = matrix(0,B,n1-n0+1)

y = matrix(rnorm(d*n),n,d)

for (u in 1:B){
  y <- y[sample(nrow(y)),]
  y.dist = as.matrix(dist(y))
  diag(y.dist) = max(y.dist)+100
  
  An = matrix(0,n,k)
  for (i in 1:n){
    An[i,] = (sort(y.dist[i,], index.return=T)$ix)[1:k]
  }


  for (t in n0:n1){
  M[u,t-n0+1] = (length(which(An[1:t,]>t))+length(which(An[(t+1):n,]<=t)))^3
  M1[u,t-n0+1] = (length(which(An[1:t,]<=t)))^3
  M2[u,t-n0+1] = (length(which(An[(t+1):n,]>t)))*(length(which(An[1:t,]<=t)))^2
  M3[u,t-n0+1] = (length(which(An[1:t,]<=t)))*(length(which(An[(t+1):n,]>t)))^2
  M4[u,t-n0+1] = (length(which(An[(t+1):n,]>t)))^3

  }
}
m = colMeans(M)
m1 = colMeans(M1)
m2 = colMeans(M2)
m3 = colMeans(M3)
m4 = colMeans(M4)
E03diff = (E3-m)/m
plot(E03diff)
E13diff = (E31-m1)/m1
plot(E13diff)
E1221diff = (E32-m2)/m2
plot(E1221diff)
E1122diff = (E33-m3)/m3
E23diff = (E34-m4)/m4
