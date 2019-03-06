k = 3
d = 100
b = 3.05
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
count = 0
double = rep(0,n)
for (i in 1:n){
  ids = An[i,]
  double[i] = length(which(An[ids,]==i))
  count = count + length(which(An[ids,]==i))
}


R2 = (cn + k)*n*k
c_t = rep(0,n1-n0+1)
E3 = rep(0,n1-n0+1)
E3hao = rep(0,n1-n0+1)
E = rep(0,n1-n0+1)
r1 = rep(0,n1-n0+1)
r2 = rep(0,n1-n0+1)
r3 = rep(0,n1-n0+1)
r4 = rep(0,n1-n0+1)
Gammahao = rep(0,n1-n0+1)
Thetahao = rep(0,n1-n0+1)
Shao = rep(0,n1-n0+1)
Var = rep(0,n1-n0+1)
Gamma = rep(0,n1-n0+1)
Theta = rep(0,n1-n0+1)
S = rep(0,n1-n0+1)
N = rep(0,23)
P = matrix(0,23,n1-n0+1)
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
  c_t[t-n0+1] = -((n - 1)*(2*R2*n + 4*count*n - R2*n^2 + R2*n^3 + 4*R2*t^2 - 4*count*n^2 + 8*count*t^2 + 4*k*n^2 - 4*k*n^3 + 6*k^2*n^2 + k^2*n^3 - k^2*n^4 - 8*count*n*t^2 + 8*count*n^2*t + 8*k*n*t^2 - 8*k*n^2*t + 8*k*n^3*t - 8*k*n^2*t^2 + 12*k^2*n*t^2 - 12*k^2*n^2*t + 4*k^2*n^3*t - 4*R2*n*t - 8*count*n*t - 4*k^2*n^2*t^2 + 4*R2*n*t^2 - 4*R2*n^2*t))/(2*t*(n - t)*(2*R2 + 4*count - 3*R2*n - 8*count*n + 4*k*n + 2*R2*n^2 - R2*n^3 + 4*R2*t^2 + 4*count*n^2 - 4*count*t^2 - 8*k*n^2 + 6*k^2*n + 4*k*n^3 - 5*k^2*n^2 - 2*k^2*n^3 + k^2*n^4 + 4*count*n*t^2 - 4*count*n^2*t - 4*k*n*t^2 + 4*k*n^2*t - 4*k*n^3*t + 4*k*n^2*t^2 - 12*k^2*n*t^2 + 12*k^2*n^2*t - 4*k^2*n^3*t - 4*R2*n*t + 4*count*n*t + 4*k^2*n^2*t^2 - 4*R2*n*t^2 + 4*R2*n^2*t))
  P[1,t-n0+1] = P[2,t-n0+1] = 2*t*(n-t)/n/(n-1)
  for (i in 3:7){
    P[i,t-n0+1] = t*(n-t)/n/(n-1)
  }
  P[8,t-n0+1] = P[9,t-n0+1] = 0
  P[10,t-n0+1] = P[11,t-n0+1] = 4*t*(n-t)*(t-1)*((n-t)-1)/n/(n-1)/(n-2)/(n-3)
  for (i in 12:15){
    P[i,t-n0+1] = t*(n-t)*((n-t-1)*(n-t-2) + (t-1)*(t-2))/n/(n-1)/(n-2)/(n-3)
  }
  for (i in 16:22){
    P[i,t-n0+1] = 2*t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  }
  P[23,t-n0+1] = 8*t*(t-1)*(t-2)*(n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  E3[t-n0+1] = sum(N*P[,t-n0+1])
  
  r1[t-n0+1] = 2*t*(n-t)/n/(n-1)
  r2[t-n0+1] = 4*t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  r3[t-n0+1] = t*(n-t)*((t-1)*(t-2) + (n-t-1)*(n-t-2))/n/(n-1)/(n-2)/(n-3)
  r4[t-n0+1] = 8*t*(t-1)*(t-2)*(n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)/(n-3)/(n-4)/(n-5)
  
  E3hao[t-n0+1] = 8*k^3*n^3*r4[t-n0+1] + 12*k^2*n^2*(r2[t-n0+1] + 3*k*(r2[t-n0+1] - 2*r4[t-n0+1])) + 
    4*k*n*(3*r2[t-n0+1] - r1[t-n0+1] + 2*r3[t-n0+1] - 4*r4[t-n0+1] + 3*k*(3*r1[t-n0+1] - 2*r2[t-n0+1] - 4*r3[t-n0+1] - 4*r4[t-n0+1]) + 8*k^2*(r3[t-n0+1] - 3*r2[t-n0+1] + 5*r4[t-n0+1])) + 
    24*count/n*(k*n^2*r4[t-n0+1] + k*n*(r1[t-n0+1] + r2[t-n0+1] - 2*r3[t-n0+1] - 4*r4[t-n0+1]) + 2*n*(2*r3[t-n0+1] - r1[t-n0+1] + 2*r4[t-n0+1])) + 
    12*(R2/n - k)*(k*n^2*(r2[t-n0+1] - 2*r4[t-n0+1]) + k*n*(2*r3[t-n0+1] - 5*r2[t-n0+1] + 8*r4[t-n0+1]) + n*(r1[t-n0+1] + r2[t-n0+1]- 2*r3[t-n0+1]-4*r4[t-n0+1])) + 
    4*(2*r3[t-n0+1] - 3*r2[t-n0+1] + 4*r4[t-n0+1])*sum(deg^3) + 
    24*(r1[t-n0+1] + r2[t-n0+1] - 2*r3[t-n0+1]-4*r4[t-n0+1])*(sum(deg*double)) + 
    24*(2*r4[t-n0+1] - r2[t-n0+1])*Mix  - 
    16*r4[t-n0+1]*(0.5*N[9] + 0.5*N[8])
  E3hao[t-n0+1] = E3hao[t-n0+1]/8
  
  
  
  E[t-n0+1] = 2*k*t*(n-t)/(n-1)
  Var[t-n0+1] = (2*t*(n - t)*(count + k*n))/(n*(n - 1)) - (4*k^2*t^2*(n - t)^2)/(n - 1)^2 + (t*(n - t)*(3*n*k^2 - 2*n*k + R2 - 2*count))/(n*(n - 1)) - (4*t*(n - t)*(t - 1)*(t - n + 1)*(k^2*n^2 - 3*k^2*n + k*n - R2 + count))/(n*(n - 1)*(n - 2)*(n - 3))
  Gamma[t-n0+1] = (E[t-n0+1]^3 + 3*E[t-n0+1]*Var[t-n0+1] - E3[t-n0+1])/((Var[t-n0+1])^1.5)
  Gammahao[t-n0+1] = (E[t-n0+1]^3 + 3*E[t-n0+1]*Var[t-n0+1] - E3hao[t-n0+1])/((Var[t-n0+1])^1.5)
  
  Theta[t-n0+1] = (-1 + (1+2*Gamma[t-n0+1]*b)^0.5)/Gamma[t-n0+1]
  S[t-n0+1] = exp(0.5*(b-Theta[t-n0+1])^2 + 1/6*Gamma[t-n0+1]*Theta[t-n0+1]^3)/((1+Gamma[t-n0+1]*Theta[t-n0+1])^0.5)
  
  Thetahao[t-n0+1] = (-1 + (1+2*Gammahao[t-n0+1]*b)^0.5)/Gammahao[t-n0+1]
  Shao[t-n0+1] = exp(0.5*(b-Thetahao[t-n0+1])^2 + 1/6*Gammahao[t-n0+1]*Thetahao[t-n0+1]^3)/((1+Gammahao[t-n0+1]*Thetahao[t-n0+1])^0.5)
  
  }
nu = function(x){
  2/x*(pnorm(x/2)-0.5)/((x/2)*pnorm(x/2)+dnorm(x/2))
}
temp = dnorm(b)*b*sum(c_t*S*nu(sqrt(2*c_t*b^2)))
integrandhao = c_t*Shao*nu(sqrt(2*c_t*b^2))
for (i in (n/2):1){
  if (is.na(integrandhao[i])) integrandhao[i] = integrandhao[i+1] - (integrandhao[i+2]-integrandhao[i+1])
  
  if (integrandhao[i]<0) integrandhao[i] = 0
}
for (i in (n/2):n){
  if (is.na(integrandhao[i])) integrandhao[i] = integrandhao[i-1] - (integrandhao[i-2]-integrandhao[i-1])
  
  if (integrandhao[i]<0) integrandhao[i] = 0
}

integrand = c_t*S*nu(sqrt(2*c_t*b^2))
for (i in (n/2):1){
  if (is.na(integrand[i])) integrand[i] = integrand[i+1] - (integrand[i+2]-integrand[i+1])
  
  if (integrand[i]<0) integrand[i] = 0
}
for (i in (n/2):n){
  if (is.na(integrand[i])) integrand[i] = integrand[i-1] - (integrand[i-2]-integrand[i-1])
  
  if (integrand[i]<0) integrand[i] = 0
}

temp = dnorm(b)*b*sum(integrand)
temphao = dnorm(b)*b*sum(integrandhao)
plot(Gamma)
plot(S)
