import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from straighten import *
import random
random.seed(1)
def rnd_case(r, q):
    s=[random.randint(1,3) for _ in range(r)]
    t=[]
    for _ in range(r):
        while True:
            x=Fr(random.randint(-3*q,3*q),q)
            if not (x.denominator==1 and x>0): break
        t.append(x)
    z=[complex(random.uniform(-.6,.6),random.uniform(-.6,.6)) for _ in range(r)]
    return s,t,z
worst=0
for trial in range(60):
    r=random.choice([1,2,3]); q=random.choice([1,2,3,4])
    s,t,z=rnd_case(r,q)
    R=NumRing(z,q)
    out=straighten_H(R,s,t)
    lhs=H_num(z,s,t); rhs=eval_out(R,out)
    err=abs(lhs-rhs)/max(1,abs(lhs)); worst=max(worst,err)
    dmax=max(len(h) for h,_ in out); wmax=max(sum(h) for h,_ in out)
    assert dmax<=r and wmax<=sum(s)
    if err>1e-9: print('FAIL',s,t,z,lhs,rhs)
print('worst rel err',worst)
# example 7.1
z=0.3; R=NumRing([z],2); out=straighten_H(R,[2],[Fr(3,2)]); print(out)
