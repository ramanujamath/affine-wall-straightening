import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from formulas import *
import random; random.seed(11)
def dangerous(R,out,canon):
    acc=defaultdict(R.zero)
    for (h,args),c in out.items():
        if h and h[0]==1:
            e0=R.arg_exp(*args[0])
            if e0==0:
                key=(h,tuple(R.arg_exp(*x) for x in args)) if canon else (h,args)
                acc[key]=R.add(acc[key],c)
    return sum(1 for v in acc.values() if not R.iszero(v))
stats=defaultdict(int)
for trial in range(300):
    r=random.choice([1,2,3,4]); q=random.choice([1,2,3,4]); N=random.choice([1,2,3,4])
    s=[random.randint(1,3) for _ in range(r)]
    if random.random()<.5: s[0]=1
    t=[]
    for _ in range(r):
        while True:
            x=Fr(random.randint(-3*q,3*q),q)
            if not(x.denominator==1 and x>0): break
        t.append(x)
    a=[random.randint(0,N-1) for _ in range(r)]
    R=CycRing(a,N,q); out=straighten_H(R,s,t)
    cls='(1,1)' if (s[0]==1 and a[0]==0) else 'allowed'
    nk=dangerous(R,out,False); nc=dangerous(R,out,True)
    stats[(cls,'cases')]+=1
    stats[(cls,'keyed_nonzero')]+= (nk>0); stats[(cls,'canon_nonzero')]+=(nc>0)
for k in sorted(stats): print(k,stats[k])
