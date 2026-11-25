import sys; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from formulas import *
R=CycRing([0,0,0],1,2)
out=straighten_H(R,[2,1,1],[0,Fr(1,2),-1])
acc=defaultdict(Fr)
for (h,args),c in out.items():
    key=(h,tuple(R.arg_exp(*x) for x in args))
    for e,v in c.items():   # M=2 group ring -> numbers: zeta^e = (-1)^e
        acc[key]+= v*(-1)**e
acc={k:v for k,v in acc.items() if v!=0}
def name(k):
    h,ex=k
    if not h: return '1'
    return 'Li_{%s}(%s)'%(','.join(map(str,h)),','.join('-1' if e else '1' for e in ex))
for k in sorted(acc,key=lambda k:(-len(k[0]),k)): print('%12s  %s'%(acc[k],name(k)))
