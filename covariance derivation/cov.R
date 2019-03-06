d = 5;
n = 1000;
k = 3;
s = 100;
t = 900;
number = 100;
VX = rep(0,number)
Co = rep(0,number)
Xs = 2*k*s*(n-s)/(n-1);
xs = rep(0,number);
Xt = 2*k*t*(n-t)/(n-1);
xt = rep(0,number);
x1t = rep(0,number);
x2t = rep(0,number);
x1s = rep(0,number);
x2s = rep(0,number);
for(l in 1:number){
  y = matrix(rnorm(d*n),n,d)
  y.dist = as.matrix(dist(y))
  diag(y.dist) = max(y.dist)+100
  An = matrix(0,n,k)
  for (i in 1:n){
    An[i,] = (sort(y.dist[i,], index.return=T)$ix)[1:k]
  }
  
  
  for(i in 1:s){
    for(j in 1:k){
      if(An[i,j]>s)
        xs[l]= xs[l]+1
    }
  }
  for(i in (s+1):n){
    for(j in 1:k)
      if(An[i,j]<=s)
        xs[l] = xs[l]+1
  }
  
  
  for(i in 1:t){
    for(j in 1:k){
      if(An[i,j]>t)
        xt[l]= xt[l]+1
    }
  }
  for(i in (t+1):n){
    for(j in 1:k)
      if(An[i,j]<=t)
        xt[l] = xt[l]+1
  }
  
  for(i in 1:t){
    for(j in 1:k){
      if(An[i,j]<=t)
        x1t[l]= x1t[l]+1
    }
  }
  
  for(i in 1:s){
    for(j in 1:k){
      if(An[i,j]<=s)
        x1s[l]= x1s[l]+1
    }
  }

  for(i in (t+1):n){
    for(j in 1:k)
      if(An[i,j]>t)
        x2t[l] = x2t[l]+1
  }
  for(i in (s+1):n){
    for(j in 1:k)
      if(An[i,j]>s)
        x2s[l] = x2s[l]+1
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
  pa = 2*s*(n-t)/n/(n-1);
  pb = s*(n-t)*(n+2*t-2*s-2)/n/(n-1)/(n-2);
  pc = 4*s*(n-t)*((n-s-2)*(t-s)+(s-1)*(n-s-1))/n/(n-1)/(n-2)/(n-3);
  
  
  VX[l] = na*pa+nb*pb+nc*pc-Xs*Xt;
  VX
  
}
VX111 = na*pa+nb*pb+nc*pc-Xs*Xt;
VX222 = (n*k+count)*(2*s*(n-t)/n/(n-1))+(3*n*k^2-2*n*k+sum(r^2)-2*count)*(s*(n-t)*(n+2*t-2*s-2)/n/(n-1)/(n-2))+(n^2*k^2-nb-na)*(4*s*(n-t)*((t-2)*(t-s)+(t-1)*(n-t-1))/n/(n-1)/(n-2)/(n-3))-(2*k*s*(n-s)/(n-1))*(2*k*t*(n-t)/(n-1));

((n*k+count)*(2*s*(n-t)/n/(n-1))+(3*n*k^2-2*n*k+sum(r^2)-2*count)*(s*(n-t)*(n+2*t-2*s-2)/n/(n-1)/(n-2))+(n^2*k^2-nb-na)*(4*s*(n-t)*((t-2)*(t-s)+(t-1)*(n-t-1))/n/(n-1)/(n-2)/(n-3))-(2*k*s*(n-s)/(n-1))*(2*k*t*(n-t)/(n-1)))


Cooov= cov(xs,xt);
a = (xs-Xs)*(xt-Xt);
Co2 = sum(a)/number;
Co1 = mean(VX)

 r1tr1s = mean(x1s*x1t)
E1t1s = na*t*(t-1)/n/(n-1)+nb*s*(s-1)*(t-2)/n/(n-1)/(n-2)+nc*s*(s-1)*(t-2)*(t-3)/n/(n-1)/(n-2)/(n-3)
r1tr2s = mean(x1t*x2s)
E1t2s = na*(t-s)*(t-s-1)/n/(n-1) + nb*((t-s)*(s-1)*(n-s-1)+(t-s)*(t-s-1)*(n-s-2))/n/(n-1)/(n-2) + nc*((s*(s-1)*(n-s)*(n-s-1)+2*s*(t-s)*(n-s-1)*(n-s-2)+(t-s)*(t-s-1)*(n-s-2)*(n-s-3))/n/(n-1)/(n-2)/(n-3))
r2t1s = mean(x2t*x1s)
E2t1s = nc*s*(s-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
r2t2s = mean(x2t*x2s)
E2t2s = na*(n-t)*(n-t-1)/n/(n-1)+nb*(n-t)*(n-t-1)*(n-s-2)/n/(n-1)/(n-2)+nc*(n-t)*(n-t-1)*(n-s-2)*(n-s-3)/n/(n-1)/(n-2)/(n-3)
