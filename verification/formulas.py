from straighten import *

def Phi(R, factors, I, pre):
    """Fourier symbol of a letter: list of (coef, h, (I,l))."""
    q=R.q; res=[]
    for (p,h),C in partial_fractions(factors).items():
        c=int(p*q)%q; kap=int(Fr(c,q)-p)
        base=R.mul(R.scal(pre,C*Fr(q)**(h-1)), R.mul(R.zpow(I,-kap),R.apow(I,c)))
        for l in range(q):
            res.append((R.mul(base,R.opow(l*c)),h,(I,l)))
    return res

def tensor(R, phis, pre, out):
    for combo in itertools.product(*phis):
        c=pre; hs=[]; args=[]
        for (cc,h,a) in combo:
            c=R.mul(c,cc); hs.append(h); args.append(a)
        key=(tuple(hs),tuple(args)); out[key]=R.add(out[key],c)

def depth_rm1(R, s, t):
    r=len(s); t=[Fr(x) for x in t]; out=defaultdict(R.zero)
    single=[Phi(R,[(t[j],s[j])],(j,j),R.one()) for j in range(r)]
    for i in range(r-1):
        tau=t[i]-t[i+1]
        for g,eps in defect(0,floor(tau)):
            if g>=0: fac=[(t[i]-g,s[i]),(t[i+1],s[i+1])]; pre=R.zpow((i,i),g)
            else:    fac=[(t[i],s[i]),(t[i+1]+g,s[i+1])]; pre=R.zpow((i+1,i+1),-g)
            M=Phi(R,fac,(i,i+1),R.one())
            tensor(R, single[:i]+[M]+single[i+2:], R.scal(pre,Fr(eps)), out)
    for a,eps in defect(0,floor(t[r-1])):
        c=R.scal(R.zpow((r-1,r-1),a), Fr(eps)/(Fr(a)-t[r-1])**s[r-1])
        tensor(R, single[:r-1], c, out)
    return out

def top_weight(R, s, t):
    r=len(s); t=[Fr(x) for x in t]; out=defaultdict(R.zero)
    Iplus=[i for i in range(r-1) if (t[i]-t[i+1]).denominator==1 and t[i]-t[i+1]>=1]
    for k in range(len(Iplus)+1):
        for B in itertools.combinations(Iplus,k):
            blocks=[]; cur=[0]
            for i in range(r-1):
                if i in B: cur.append(i+1)
                else: blocks.append(cur); cur=[i+1]
            blocks.append(cur)
            pre=R.one(); phis=[]
            for b in blocks:
                last=b[-1]
                for l in b: pre=R.mul(pre,R.zpow((l,l),int(t[l]-t[last])))
                phis.append(Phi(R,[(t[last],sum(s[l] for l in b))],(b[0],last),R.one()))
            tensor(R,phis,pre,out)
    return out

def part(out, pred):
    return {k:v for k,v in out.items() if pred(k)}

def same(R, A, B):
    keys=set(A)|set(B); bad=0
    for k in keys:
        d=R.add(A.get(k,R.zero()), R.scal(B.get(k,R.zero()),Fr(-1)))
        if not R.iszero(d): bad+=1
    return bad
