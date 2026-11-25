import sys,time; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from formulas import *
cases=[([2,1,1],[0,Fr(1,2),-1],2),([2,1,1],[Fr(1,3),Fr(-5,3),Fr(7,3)],3),
       ([3,2,1,1],[Fr(1,2),Fr(-3,2),Fr(5,2),-2],2),([2,2,2],[Fr(7,4),Fr(-9,4),Fr(5,4)],4),
       ([2,1,1,1],[0,-1,-2,-3],1)]
for s,t,q in cases:
    r=len(s); R=CycRing([0]*r,1,q); STATS.clear(); t0=time.time()
    out=straighten_H(R,s,t); dt=time.time()-t0
    canon=set((h,tuple(R.arg_exp(*a) for a in args)) for (h,args),c in out.items())
    w=sum(s); bound=sum(math.comb(w,d)*(q*r*(r+1)//2)**d for d in range(r+1))
    first=1
    for j in range(r): first*=1+abs(math.floor(t[j]))
    for i in range(r-1): first*=1+abs(math.floor(t[i]-t[i+1]))
    print(s,[str(x) for x in t],'q=%d'%q,'symbols=%d'%len(out),'distinct=%d'%len(canon),'fourier=%d'%STATS['fourier'],'first-layer<=%d'%first,'bound=%d'%bound,'%.2fs'%dt)
# staircase S_4: top-weight part should be zeta*(2,1,1,1) symbolically
R=CycRing([0]*4,1,1); out=straighten_H(R,[2,1,1,1],[0,-1,-2,-3])
print({ (h): c for (h,a),c in out.items() if sum(h)==5})
