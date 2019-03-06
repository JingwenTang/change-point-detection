syms k
syms t
syms n
syms s
syms count
syms R2
X0 = 2*k*t*(n-t)/(n-1);
na = n*k+count;
nb = 3*n*k^2-2*n*k+R2-2*count;
nc = n^2*k^2-nb-na;
pa = 2*t*(n-t)/(n-1)/n;
pb = t*(n-t)/n/(n-1);
pc = 4*t*(n-t)*(t-1)*((n-t)-1)/n/(n-1)/(n-2)/(n-3);
VX = na*pa+nb*pb+nc*pc-X0^2
simplify(VX)