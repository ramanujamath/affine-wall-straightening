import numpy as np, math
G=0.915965594177219015054603514932
z2=math.pi**2/6; z3=1.2020569031595942854; z4=math.pi**4/90; z5=1.0369277551433699263
def S(f,N): n=np.arange(1,N+1,dtype=np.float64); return np.sum(f(n)[::-1])
def ext(f):  # Richardson assuming tail ~ a/N + b/N^2 + c/N^3
    Ns=[200000,400000,800000,1600000]; v=[S(f,N) for N in Ns]
    for k in (1,2,3):
        v=[(2**k*v[i+1]-v[i])/(2**k-1) for i in range(len(v)-1)]
    return v[0]
# 7.1
z=0.3;a=math.sqrt(z)
Li2=lambda x: sum(x**k/k**2 for k in range(1,400))
print('7.1', sum(z**n/(n-1.5)**2 for n in range(1,400)) - (4*z+2*z*a*(Li2(a)-Li2(-a))))
r3=math.sqrt(3);l3=math.log(3);pi=math.pi
print('7.3',ext(lambda n:1/(n*(n-2/3)))-(3*pi/(4*r3)+9/4*l3))
print('7.4',ext(lambda n:1/(n**2*(n-2/3)))-(-1.5*z2+9*pi/(8*r3)+27/8*l3))
print('7.5',ext(lambda n:1/((n-.25)**2*(n-.75)))-(4*pi-2*pi**2+16*G))
print('7.6',ext(lambda n:1/(n**2*(n-1/3)*(n-2/3)))-(4.5*z2+135*pi/(8*r3)-243/8*l3))
