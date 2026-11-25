import numpy as np, math, itertools
N=3000
def excl_cumsum(v): c=np.cumsum(v); return np.concatenate([[0.0],c[:-1]])  # sum over indices < n
n=np.arange(1,N+1,dtype=float)
# X with outer n1<=N
i3=excl_cumsum(1/(n+1)); i2=excl_cumsum(i3/(n-0.5)); X=np.sum(i2/n**2)
# F: Li_{2,1,1}(x,y,w) truncated n1<=2N
M=2*N; k=np.arange(1,M+1,dtype=float)
def Li211(x,y,w):
    a=excl_cumsum(w**k/k); b=excl_cumsum(a*y**k/k); return np.sum(b*x**k/k**2)
F=2*sum((-1)**b*Li211((-1.)**a,(-1.)**b,(-1.)**c) for a,b,c in itertools.product((0,1),repeat=3))
# C12: m>n>0, 1/(m^2 (m-1/2) n)
C12=np.sum(excl_cumsum(1/n)/(n**2*(n-0.5)))
# W3: m1>m2>1, 1/(m1^2 (m2-1/2))
v=1/(n-0.5); v[0]=0; W3=np.sum(excl_cumsum(v)/n**2)
W23=np.sum(excl_cumsum(1/((n-0.5)*n))/n**2)
W323=2*np.sum(1/n[1:]**2)
print('X',X,'F',F,'C12',C12,'W3',W3,'W23',W23); print('resid', X-(F-C12-W3+W23-W323))
# staircase S3,S4 via nested sums with Richardson in N
z2=math.pi**2/6; z3=1.2020569031595942854; z4=math.pi**4/90; z5=1.0369277551433699263
def Sr(r,N):
    n=np.arange(1,N+1,dtype=float); inner=np.ones(N)
    for j in range(r,1,-1): inner=excl_cumsum(inner/(n+j-1))
    return np.sum(inner/n**2)
for r,val in [(3,3*z4-3*z3+z2/2),(4,4*z5-5.5*z4+2*z3-z2/6)]:
    print('S%d'%r,[Sr(r,N)-val for N in (10**5,10**6,4*10**6)])
