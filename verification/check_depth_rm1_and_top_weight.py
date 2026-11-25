import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from formulas import *
import random; random.seed(7)
nd=nw=0
for trial in range(80):
    r=random.choice([2,3,4]); q=random.choice([1,2,3])
    s=[random.randint(1,3) for _ in range(r)]
    t=[]
    for _ in range(r):
        while True:
            x=Fr(random.randint(-3*q,3*q),q)
            if not(x.denominator==1 and x>0): break
        t.append(x)
    a=[random.randint(0,5) for _ in range(r)]; N=6
    R=CycRing(a,N,q)            # exact arithmetic, roots of unity
    out=straighten_H(R,s,t)
    w=sum(s)
    A=part(out,lambda k:len(k[0])==r-1); B=dict(depth_rm1(R,s,t))
    nd+=same(R,A,B)
    A=part(out,lambda k:sum(k[0])==w); B=dict(top_weight(R,s,t))
    nw+=same(R,A,B)
print('depth r-1 mismatches:',nd,' top-weight mismatches:',nw)
