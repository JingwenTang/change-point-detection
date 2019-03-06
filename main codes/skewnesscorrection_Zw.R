k = 3
d = 100
b = 3.6
n = 1000
n0 = 25
n1 = n-n0

y = matrix(rnorm(d*n),n,d)
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
mix = 0
Mix = 0
for (i in 1:n){
  for (j in 1:k){
    Mix = Mix + deg[i]*deg[An[i,j]]
    mix = mix + deg[i]*(deg[An[i,j]]-1)
  }
}
cn = sum((deg-k)^2)/n/k
R2 = (cn + k)*n*k
count = 0
double = rep(0,n)
for (i in 1:n){
  ids = An[i,]
  double[i] = length(which(An[ids,]==i))
  count = count + length(which(An[ids,]==i))
}
R2 = (cn + k)*n*k
ct = rep(0,n1-n0+1)
E3w = rep(0,n1-n0+1)
q = rep(0,n1-n0+1)
p = rep(0,n1-n0+1)
E = rep(0,n1-n0+1)
Var = rep(0,n1-n0+1)
Gamma = rep(0,n1-n0+1)
Theta = rep(0,n1-n0+1)
K = rep(0,n1-n0+1)
E31 = rep(0,n1-n0+1)
E32 = rep(0,n1-n0+1)
E33 = rep(0,n1-n0+1)
E34 = rep(0,n1-n0+1)
P1 = matrix(0,23,n1-n0+1)
P2 = matrix(0,23,n1-n0+1)
P3 = matrix(0,23,n1-n0+1)
P4 = matrix(0,23,n1-n0+1)
N = rep(0,23)
N[1] = n*k
N[2] = 3*count
N[3] = 6*n*k^2 - 6*count
N[4] = 3*n*k*(k-1)
N[5] = 3*R2 - 3*n*k
N[6] = 6*(k-1)*count
N[7] = 6*(sum(deg*double)-count)
for (i in 1:n){
  for (j in 1:k){
    for (m in 1:k){
      N[8] = N[8] + length(which(An[An[i,j],]==An[i,m]))
      N[9] = N[9] + length(which(An[An[An[i,j],m],]==i))
    }
  }
}
N[8] = 6*N[8]
N[9] = 2*N[9]
N[10] = 3*(n^2*k^2 - 3*n*k^2 + n*k + count - R2)
N[11] = 3*count*(n*k - 2*k) - N[7]
N[12] = n*k*(k-1)*(k-2)
N[13] = 3*n*k^2*(k-1) - N[6]
N[14] = 3*k*R2 - 3*n*k^2 - N[7]
N[15] = sum(deg^3) - 3*R2 + 2*n*k
N[16] = 6*n*k^3 - N[6] - N[7] - 2*N[2] - 3*N[9]
N[17] = 6*mix - N[7] - N[8]
N[18] = 6*n*k^2*(k-1) - N[6] - N[8]
N[19] = 6*(k-1)*(R2 - n*k) - N[8]
N[20] = 3*k^2*(k-1)*(n-3)*n - N[13] - N[19]
N[21] = 3*k*(n-3)*(R2 - n*k) - 3*N[15] - N[17]
N[22] = 6*(n-3)*k*(n*k^2 - count) - 2*N[14] - N[16] - N[17]
N[23] = n^3*k^3 - sum(N[1:22])


for (t in n0:n1){
  ct[t-n0+1] = -((n - 1)*(2*t^2 - 2*n*t + n))/(2*t*(t - 1)*(n^2 - 2*n*t - n + t^2 + t))
  
  P1[1,t-n0+1] = P1[2,t-n0+1] = t*(t-1)/n/(n-1)
  for (i in 3:9){
    P1[i,t-n0+1] = t*(t-1)*(t-2)/n/(n-1)/(n-2)
  }
  
  for (i in 10:19){
    P1[i,t-n0+1] = t*(t-1)*(t-2)*(t-3)/n/(n-1)/(n-2)/(n-3)
  }
  for (i in 20:22){
    P1[i,t-n0+1] = t*(t-1)*(t-2)*(t-3)*(t-4)/n/(n-1)/(n-2)/(n-3)/(n-4)
  }
  P1[23,t-n0+1] = t*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  
  for (i in 1:9){
    P2[i,t-n0+1] = 0
  }
  for (i in 10:11){
    P2[i,t-n0+1] = 1/3*t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  }
  for (i in 12:19){
    P2[i,t-n0+1] = 0
  }
  for (i in 20:22){
    P2[i,t-n0+1] = 1/3*t*(t-1)*(t-2)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)/(n-4)
  }
  P2[23,t-n0+1] = t*(t-1)*(t-2)*(t-3)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  
  for (i in 1:9){
    P3[i,t-n0+1] = 0
  }
  for (i in 10:11){
    P3[i,t-n0+1] = 1/3*t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  }
  for (i in 12:19){
    P3[i,t-n0+1] = 0
  }
  for (i in 20:22){
    P3[i,t-n0+1] = 1/3*t*(t-1)*(n-t-2)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)/(n-4)
  }
  P3[23,t-n0+1] = t*(t-1)*(n-t-2)*(n-t-3)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  
  P4[1,t-n0+1] = P4[2,t-n0+1] = (n-t)*((n-t)-1)/n/(n-1)
  for (i in 3:9){
    P4[i,t-n0+1] = (n-t)*((n-t)-1)*((n-t)-2)/n/(n-1)/(n-2)
  }
  
  for (i in 10:19){
    P4[i,t-n0+1] = (n-t)*((n-t)-1)*((n-t)-2)*((n-t)-3)/n/(n-1)/(n-2)/(n-3)
  }
  for (i in 20:22){
    P4[i,t-n0+1] = (n-t)*((n-t)-1)*((n-t)-2)*((n-t)-3)*((n-t)-4)/n/(n-1)/(n-2)/(n-3)/(n-4)
  }
  P4[23,t-n0+1] = (n-t)*((n-t)-1)*((n-t)-2)*((n-t)-3)*((n-t)-4)*((n-t)-5)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  
  
  E31[t-n0+1] = sum(N*P1[,t-n0+1])
  E32[t-n0+1] = sum(N*P2[,t-n0+1])
  E33[t-n0+1] = sum(N*P3[,t-n0+1])
  E34[t-n0+1] = sum(N*P4[,t-n0+1])
  q[t-n0+1] = (n-t-1)/(n-2)
  p[t-n0+1] = (t-1)/(n-2)
  E3w[t-n0+1] = (q[t-n0+1])^3*E31[t-n0+1] + 3*(q[t-n0+1])^2*p[t-n0+1]*E32[t-n0+1] + 3*q[t-n0+1]*(p[t-n0+1])^2*E33[t-n0+1] + (p[t-n0+1])^3*E34[t-n0+1]
  
  E[t-n0+1] = -(k*n*(t - 1)*(t - n + 1))/(n^2 - 3*n + 2)
  Var[t-n0+1] = (t*(t - 1)*(n^2 - 2*n*t - n + t^2 + t)*(R2 + 2*count - R2*n - 3*count*n + 2*k*n + count*n^2 - 3*k*n^2 + 3*k^2*n + k*n^3 - k^2*n^2))/(n*(n - 3)*(n^2 - 3*n + 2)^2)
  Gamma[t-n0+1] = (E3w[t-n0+1] - 3*E[t-n0+1]*Var[t-n0+1] - E[t-n0+1]^3)/((Var[t-n0+1])^1.5)
  Theta[t-n0+1] = (-1 + (1+2*Gamma[t-n0+1]*b)^0.5)/Gamma[t-n0+1]
  K[t-n0+1] = exp(0.5*(b-Theta[t-n0+1])^2 + 1/6*Gamma[t-n0+1]*Theta[t-n0+1]^3)/((1+Gamma[t-n0+1]*Theta[t-n0+1])^0.5)
  
}               
nu = function(x){
  2/x*(pnorm(x/2)-0.5)/((x/2)*pnorm(x/2)+dnorm(x/2))
}
integrand = ct*K*nu(sqrt(2*ct*b^2))
for (i in floor((n1-n0+1)/2):1){
  if (is.na(integrand[i])) integrand[i] = integrand[i+1] - (integrand[i+2]-integrand[i+1])
  
  if (integrand[i]<0) integrand[i] = 0
}
for (i in ceiling((n1-n0+1)/2):(n1-n0+1)){
  if (is.na(integrand[i])) integrand[i] = integrand[i-1] - (integrand[i-2]-integrand[i-1])
  
  if (integrand[i]<0) integrand[i] = 0
}

temp = dnorm(b)*b*sum(integrand)
plot(K)
plot(Gamma)

