d = 5;
n = 1000;
k = 3;
n1 = 199;
n2 = n-n1;
number = 100;
VX = rep(0,number)
X0 = 2*k*n1*n2/(n-1)
X1 = k*n1*(n1-1)/(n-1)
X2 = k*n2*(n2-1)/(n-1)


x0 = rep(0,number)
x1 = rep(0,number)
x2 = rep(0,number)
S = rep(0,number)

for(l in 1:number){
  y = matrix(rnorm(d*n),n,d)
y.dist = as.matrix(dist(y))
diag(y.dist) = max(y.dist)+100
An = matrix(0,n,k)
for (i in 1:n){
  An[i,] = (sort(y.dist[i,], index.return=T)$ix)[1:k]
}


for(i in 1:n1){
  for(j in 1:k){
    if(An[i,j]>n1)
      x0[l]= x0[l]+1
  }
}
for(i in (n1+1):n){
  for(j in 1:k)
    if(An[i,j]<=n1)
      x0[l] = x0[l]+1
}


for(i in 1:n1){
  for(j in 1:k){
    if(An[i,j]<=n1)
      x1[l]= x1[l]+1
  }
}
for(i in (n1+1):n){
  for(j in 1:k)
    if(An[i,j]>n1)
      x2[l] = x2[l]+1
}




temp = table(An)
id = as.numeric(row.names(temp))
deg = rep(0,n)
deg[id] = temp
count = 0
for (i in 1:n){
  ids = An[i,]
  count = count + length(which(An[ids,]==i))
}


r = deg
na = n*k+count;
nb = 3*n*k^2-2*n*k+sum(r^2)-2*count;
nc = n^2*k^2-nb-na;
pa = 2*n1*n2/(n-1)/n;
pb = n1*n2/n/(n-1);
pc = 4*n1*n2*(n1-1)*(n2-1)/n/(n-1)/(n-2)/(n-3);
VX[l] = na*pa+nb*pb+nc*pc-X0^2;

t = n1
ER1R2 = nc*n1*n2*(n1-1)*(n2-1)/n/(n-1)/(n-2)/(n-3);
covr1r2 = ER1R2 - X1*X2
varianceR1 = n1/(n - 1)*(n1 - 1)*(count + k*n)/n +nb*n1*(n1-1)*(n1-2)/n/(n-1)/(n-2) + nc*n1*(n1-1)*(n1-2)*(n1-3)/n/(n-1)/(n-2)/(n-3) - X1^2
varianceR2 = n2/(n - 1)*(n2 - 1)*(count + k*n)/n +nb*n2*(n2-1)*(n2-2)/n/(n-1)/(n-2) + nc*n2*(n2-1)*(n2-2)*(n2-3)/n/(n-1)/(n-2)/(n-3) - X2^2
temp = matrix(c(varianceR1,covr1r2,covr1r2,varianceR2),nrow = 2,ncol = 2,byrow = TRUE)
M = solve(temp)
temp = c(x1[l]-X1,x2[l] - X2)%*%M
S[l] = temp%*%matrix(x1[l] - X1, x2[l] - X2,nrow = 2,ncol = 1)
}



R0 = mean(x0)
E33 = mean(x0^3)
R1 = mean(x1)
R2 = mean(x2)
V0 = var(x0)
V1 = var(x1)
V2 = var(x2)



Variance = mean(VX)

cn = sum((deg-k)^2)/n/k

vn = count/n/k
t = n1
EX = 4*k*t*(n-t)/(n-1)
h = 4*(t-1)*(n-t-1)/((n-2)*(n-3))
VX1 = EX*(h*(1+vn-2*k/(n-1))+(1-h)*cn)/4
(2*t*(n - t)*(count + k*n))/(n*(n - 1)) - (4*k^2*t^2*(n - t)^2)/(n - 1)^2 + (t*(n - t)*(3*n*k^2 - 2*n*k + R2 - 2*count))/(n*(n - 1)) - (4*t*(n - t)*(t - 1)*(t - n + 1)*(k^2*n^2 - 3*k^2*n + k*n - R2 + count))/(n*(n - 1)*(n - 2)*(n - 3))




r = deg
Vat = ((4*k*t*(n-t)/(n-1))*((4*(t-1)*(n-t-1)/((n-2)*(n-3)))*(1+(count/n/k)-2*k/(n-1))+(1-(4*(t-1)*(n-t-1)/((n-2)*(n-3))))*(sum(r^2)/k/n-k))/4)^(-0.5)



r = deg
na = n*k+count;
nb = 3*n*k^2-2*n*k+sum(r^2)-2*count;
nc = n^2*k^2-nb-na;
pa = 2*n1*n2/(n-1)/n;
pb = n1*n2/n/(n-1);
pc = 4*n1*n2*(n1-1)*(n2-1)/n/(n-1)/(n-2)/(n-3);
VX2 = na*pa+nb*pb+nc*pc-X0^2
t = n1
n2 = n-n1
R2 = sum(r^2)
(2*t*(n - t)*(count + k*n))/(n*(n - 1)) - (4*k^2*t^2*(n - t)^2)/(n - 1)^2 + (t*(n - t)*(3*n*k^2 - 2*n*k + R2 - 2*count))/(n*(n - 1)) - (4*t*(n - t)*(t - 1)*(t - n + 1)*(k^2*n^2 - 3*k^2*n + k*n - R2 + count))/(n*(n - 1)*(n - 2)*(n - 3))




varianceR1 = n1/(n - 1)*(n1 - 1)*(count + k*n)/n +nb*n1*(n1-1)*(n1-2)/n/(n-1)/(n-2) + nc*n1*(n1-1)*(n1-2)*(n1-3)/n/(n-1)/(n-2)/(n-3) - X1^2
varianceR2 = n2/(n - 1)*(n2 - 1)*(count + k*n)/n +nb*n2*(n2-1)*(n2-2)/n/(n-1)/(n-2) + nc*n2*(n2-1)*(n2-2)*(n2-3)/n/(n-1)/(n-2)/(n-3) - X2^2

(((4*k*t*(n-t)/(n-1))*((4*(t-1)*(n-t-1)/((n-2)*(n-3)))*(1+(count/n/k)-2*k/(n-1))+(1-(4*(t-1)*(n-t-1)/((n-2)*(n-3))))*(sum(r^2)/k/n-k))/4)^(-0.5))
(((4*k*s*(n-s)/(n-1))*((4*(s-1)*(n-s-1)/((n-2)*(n-3)))*(1+(count/n/k)-2*k/(n-1))+(1-(4*(s-1)*(n-s-1)/((n-2)*(n-3))))*(sum(r^2)/k/n-k))/4)^(-0.5))


t = n1
ER1R2 = nc*n1*n2*(n1-1)*(n2-1)/n/(n-1)/(n-2)/(n-3);
covr1r2 = ER1R2 - X1*X2
a1 = (n-t-1)/(n-2)
a2 = (t-1)/(n-2)
ERw = a1*X1 + a2*X2
ERdiff = X1-X2
VRw = a1^2*varianceR1 + a2^2*varianceR2 + 2*a1*a2*covr1r2
VRdiff = varianceR1 + varianceR2 - 2*covr1r2

