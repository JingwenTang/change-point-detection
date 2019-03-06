


VX = expression((n*k+count)*(2*s*(n-t)/n/(n-1))+(3*n*k^2-2*n*k+sum(r^2)-2*count)*(s*(n-t)*(n+2*t-2*s-2)/n/(n-1)/(n-2))+(n^2*k^2-nb-na)*(4*s*(n-t)*((t-2)*(t-s)+(t-1)*(n-t-1))/n/(n-1)/(n-2)/(n-3))-(2*k*s*(n-s)/(n-1))*(2*k*t*(n-t)/(n-1)));

Vat = expression(((4*k*t*(n-t)/(n-1))*((4*(t-1)*(n-t-1)/((n-2)*(n-3)))*(1+(count/n/k)-2*k/(n-1))+(1-(4*(t-1)*(n-t-1)/((n-2)*(n-3))))*(sum(r^2)/k/n-k))/4)^(-0.5))


Vas = expression(((4*k*s*(n-s)/(n-1))*((4*(s-1)*(n-s-1)/((n-2)*(n-3)))*(1+(count/n/k)-2*k/(n-1))+(1-(4*(s-1)*(n-s-1)/((n-2)*(n-3))))*(sum(r^2)/k/n-k))/4)^(-0.5))

Cov = expression((((4*k*t*(n-t)/(n-1))*((4*(t-1)*(n-t-1)/((n-2)*(n-3)))*(1+(count/n/k)-2*k/(n-1))+(1-(4*(t-1)*(n-t-1)/((n-2)*(n-3))))*(R2/k/n-k))/4)^(-0.5))*(((4*k*s*(n-s)/(n-1))*((4*(s-1)*(n-s-1)/((n-2)*(n-3)))*(1+(count/n/k)-2*k/(n-1))+(1-(4*(s-1)*(n-s-1)/((n-2)*(n-3))))*(R2/k/n-k))/4)^(-0.5))*((n*k+count)*(2*s*(n-t)/n/(n-1))+(3*n*k^2-2*n*k+R2-2*count)*(s*(n-t)*(n+2*t-2*s-2)/n/(n-1)/(n-2))+(n^2*k^2-(3*n*k^2-2*n*k+R2-2*count)-(n*k+count))*(4*s*(n-t)*((t-2)*(t-s)+(t-1)*(n-t-1))/n/(n-1)/(n-2)/(n-3))-(2*k*s*(n-s)/(n-1))*(2*k*t*(n-t)/(n-1))));

c = D(Cov,"s")


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
cn = sum((deg-k)^2)/n/k
count = 0
for (i in 1:n){
  ids = An[i,]
  count = count + length(which(An[ids,]==i))
}
R2 = (cn + k)*n*k
c_t = rep(0,n1-n0+1)
CT = rep(0,n1-n0+1)
ct = rep(0,n1-n0+1)
miu = rep(0,n1-n0+1)
MIU1 = rep(0,n1-n0+1)
sigma = rep(0,n1-n0+1)
MIU = rep(0,n1-n0+1)
MIU2 = rep(0,n1-n0+1)
decimal = rep(0,n1-n0+1)
for (t in n0:n1){



c_t[t-n0+1] = -((n - 1)*(2*R2*n + 4*count*n - R2*n^2 + R2*n^3 + 4*R2*t^2 - 4*count*n^2 + 8*count*t^2 + 4*k*n^2 - 4*k*n^3 + 6*k^2*n^2 + k^2*n^3 - k^2*n^4 - 8*count*n*t^2 + 8*count*n^2*t + 8*k*n*t^2 - 8*k*n^2*t + 8*k*n^3*t - 8*k*n^2*t^2 + 12*k^2*n*t^2 - 12*k^2*n^2*t + 4*k^2*n^3*t - 4*R2*n*t - 8*count*n*t - 4*k^2*n^2*t^2 + 4*R2*n*t^2 - 4*R2*n^2*t))/(2*t*(n - t)*(2*R2 + 4*count - 3*R2*n - 8*count*n + 4*k*n + 2*R2*n^2 - R2*n^3 + 4*R2*t^2 + 4*count*n^2 - 4*count*t^2 - 8*k*n^2 + 6*k^2*n + 4*k*n^3 - 5*k^2*n^2 - 2*k^2*n^3 + k^2*n^4 + 4*count*n*t^2 - 4*count*n^2*t - 4*k*n*t^2 + 4*k*n^2*t - 4*k*n^3*t + 4*k*n^2*t^2 - 12*k^2*n*t^2 + 12*k^2*n^2*t - 4*k^2*n^3*t - 4*R2*n*t + 4*count*n*t + 4*k^2*n^2*t^2 - 4*R2*n*t^2 + 4*R2*n^2*t))


}

nu = function(x){
  2/x*(pnorm(x/2)-0.5)/((x/2)*pnorm(x/2)+dnorm(x/2))
}
p = dnorm(b)/b*sum(sigma)

temp = dnorm(b)*b*sum(c_t*nu(sqrt(2*c_t*b^2)))


plot(n0:n1,c_t)
plot(n0:n1,ct)
plot(n0:n1,ct-c_t)
min(c_t)
 n
 t = 5000
 
 2*t*(n - t)*(2*R2 + 4*count - 3*R2*n - 8*count*n + 4*k*n + 2*R2*n^2 - R2*n^3 + 4*R2*t^2 + 4*count*n^2 - 4*count*t^2 - 8*k*n^2 + 6*k^2*n + 4*k*n^3 - 5*k^2*n^2 - 2*k^2*n^3 + k^2*n^4 + 4*count*n*t^2 - 4*count*n^2*t - 4*k*n*t^2 + 4*k*n^2*t - 4*k*n^3*t + 4*k*n^2*t^2 - 12*k^2*n*t^2 + 12*k^2*n^2*t - 4*k^2*n^3*t - 4*R2*n*t + 4*count*n*t + 4*k^2*n^2*t^2 - 4*R2*n*t^2 + 4*R2*n^2*t)
(2*R2*n + 4*count*n - R2*n^2 + R2*n^3 + 4*R2*t^2 - 4*count*n^2 + 8*count*t^2 + 4*k*n^2 - 4*k*n^3 + 6*k^2*n^2 + k^2*n^3 - k^2*n^4 - 8*count*n*t^2 + 8*count*n^2*t + 8*k*n*t^2 - 8*k*n^2*t + 8*k*n^3*t - 8*k*n^2*t^2 + 12*k^2*n*t^2 - 12*k^2*n^2*t + 4*k^2*n^3*t - 4*R2*n*t - 8*count*n*t - 4*k^2*n^2*t^2 + 4*R2*n*t^2 - 4*R2*n^2*t)
 