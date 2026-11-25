import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from collections import defaultdict
from straighten import CycRing, straighten_H
R=CycRing([0,0,0],1,2); out=straighten_H(R,[2,1,1],[0,Fr(1,2),-1]); acc=defaultdict(Fr)
for (h,args),c in out.items():
    key=(h,tuple(R.arg_exp(*x) for x in args))
    for e,v in c.items(): acc[key]+=v*(-1)**e
acc={k:v for k,v in acc.items() if v!=0}
def nested(letters,N):
    n=np.arange(1,N+2,dtype=float); inner=np.ones(N+1)
    for sgn,f in reversed(letters):
        term=(sgn**np.arange(1,N+2))*f(n)*inner
        inner=np.concatenate([[0],np.cumsum(term)[:-1]])
    c=np.cumsum(term); return 0.5*(c[N-1]+c[N])   # average of partial sums N, N+1
def rhs(N):
    tot=0.
    for (h,ex),v in acc.items():
        if not h: tot+=float(v); continue
        L=[((-1.) if e else 1., (lambda hh: (lambda n: n**-hh))(hh)) for hh,e in zip(h,ex)]
        tot+=float(v)*nested(L,N)
    return tot
def X(N):
    return nested([(1.,lambda n:n**-2.),(1.,lambda n:1/(n-.5)),(1.,lambda n:1/(n+1))],N)
for N in [10**5,10**6,4*10**6]:
    print(N, X(N), rhs(N), X(N)-rhs(N))
Ns=[2*10**5,5*10**5,10**6,2*10**6,4*10**6,8*10**6]
xs=np.array([X(N) for N in Ns]); rs=np.array([rhs(N) for N in Ns])
def fit(vals):
    A=np.array([[1,np.log(N)**2/N,np.log(N)/N,1/N] for N in Ns]); return np.linalg.lstsq(A,vals,rcond=None)[0][0]
print('X_inf  ~',fit(xs)); print('RHS_inf~',fit(rs)); print('diff   ~',fit(xs)-fit(rs))
