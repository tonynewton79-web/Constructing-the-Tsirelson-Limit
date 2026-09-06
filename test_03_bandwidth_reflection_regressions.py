#!/usr/bin/env python3
"""
GitHub test 3/3 — numerical regressions, reflection branch, and bandwidth attack.

This file consolidates:
  * independent high-precision reconstruction of the sigma=0.02 covariance;
  * L=8 Fejer norm and finite-Weyl Bell values;
  * reflection correlation and Fejer-reflection benchmark;
  * exact rational Sturm bound + numerical value for the optimized hybrid axis;
  * end-to-end exact sigma=0.01 witness certificate;
  * deterministic degree-15 bandwidth continuation to K(0)=I_8.

Run:
    python test_03_bandwidth_reflection_regressions.py

Optional pytest:
    pytest -q test_03_bandwidth_reflection_regressions.py

Dependencies:
    mpmath
    numpy
    scipy
"""
from fractions import Fraction as F
import mpmath as mp

def _trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p
def _der(p): return [F(i)*p[i] for i in range(1,len(p))] or [F(0)]
def _divrem(a,b):
    a=_trim(a);b=_trim(b);r=a[:]
    while len(r)>=len(b) and not (len(r)==1 and r[0]==0):
        c=r[-1]/b[-1];d=len(r)-len(b)
        for i in range(len(b)):r[d+i]-=c*b[i]
        r=_trim(r)
    return r
def _sturm(p):
    seq=[_trim(p),_der(_trim(p))]
    while True:
        r=_divrem(seq[-2],seq[-1])
        if len(r)==1 and r[0]==0:break
        seq.append([-x for x in r])
    return seq
def _eval(p,x):
    s=F(0)
    for c in reversed(p):s=s*x+c
    return s
def _var(vals):
    s=[1 if v>0 else -1 for v in vals if v]
    return sum(s[i]!=s[i-1] for i in range(1,len(s)))
def _padd(p,q):
    n=max(len(p),len(q));out=[F(0)]*n
    for i in range(n):out[i]=(p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0)
    return _trim(out)
def _pscale(p,s):return [s*x for x in p]
def _px(p):return [F(0)]+list(p)
def _Tpoly(n):
    if n==0:return [F(1)]
    if n==1:return [F(0),F(1)]
    a,b=[F(1)],[F(0),F(1)]
    for _ in range(1,n):a,b=b,_padd(_pscale(_px(b),2),_pscale(a,-1))
    return b

def _Tn(n,y):
    if n==0:return mp.mpf(1)
    if n==1:return y
    a,b=mp.mpf(1),y
    for _ in range(1,n):a,b=b,2*y*b-a
    return b
def _Un(n,y):
    if n==0:return mp.mpf(1)
    if n==1:return 2*y
    a,b=mp.mpf(1),2*y
    for _ in range(1,n):a,b=b,2*y*b-a
    return b

def _continuous_norm(coeff,ms):
    def P(y):return mp.fsum(coeff[i]*_Tn(ms[i],y) for i in range(len(ms)))
    def D(y):return mp.fsum(coeff[i]*ms[i]*_Un(ms[i]-1,y) for i in range(len(ms)))
    roots=[];N=20000
    x0=mp.mpf(0);f0=D(x0)
    for i in range(1,N+1):
        x1=mp.mpf(i)/N;f1=D(x1)
        if f0*f1<0:
            a,b=x0,x1;fa=f0
            for _ in range(100):
                m=(a+b)/2;fm=D(m)
                if fa*fm<=0:b=m
                else:a=m;fa=fm
            roots.append((a+b)/2)
        x0,f0=x1,f1
    vals=[abs(P(mp.mpf(0))),abs(P(mp.mpf(1)))]
    vals.extend(abs(P(r)) for r in roots)
    return max(vals),roots

def run_sigma002_regressions():
    mp.mp.dps=70
    sigma=mp.mpf("0.02")
    def coth(x):return 1/mp.tanh(x)
    A0=mp.quad(lambda k:k*coth(mp.pi*k/2)*mp.e**(-(k/sigma)**2),[0,mp.inf])
    B0=mp.quad(lambda k:k*mp.tanh(mp.pi*k/2)*mp.e**(-(k/sigma)**2),[0,mp.inf])
    Bc=mp.pi**mp.mpf("1.5")*sigma**3/8*mp.e**(mp.pi**2*sigma**2/4)
    Ac=A0+Bc
    nuA=mp.sqrt(Ac/Bc);nuB=mp.sqrt(A0/B0)
    kq=A0**mp.mpf("0.75")/(Ac*Bc*B0)**mp.mpf("0.25")
    kp=(Ac*A0*B0)**mp.mpf("0.25")/Bc**mp.mpf("0.75")
    expected=[
        (nuA,mp.mpf("45.00841733645926200048")),
        (nuB,mp.mpf("45.03061897035379099827")),
        (kq,mp.mpf("45.00840364042463202550")),
        (kp,mp.mpf("44.98621295944687856821")),
    ]
    for got,ref in expected:assert abs(got-ref)<mp.mpf("5e-20")

    L=8;ms=[2*r+1 for r in range(L)]
    coeff=[4/mp.pi*((-1)**r)/ms[r]*(1-mp.mpf(ms[r])/(2*L)) for r in range(L)]
    M8=mp.fsum(coeff) # monotone Fejer polynomial on [0,1]
    assert abs(M8-mp.mpf("0.96036378670045260128"))<mp.mpf("5e-20")

    def corr(c,kappa):
        s=mp.mpf(0)
        for r,m in enumerate(ms):
            for t,n in enumerate(ms):
                em=mp.e**(-mp.pi/4*(m*m*nuA+n*n*nuB-2*m*n*kappa))
                ep=mp.e**(-mp.pi/4*(m*m*nuA+n*n*nuB+2*m*n*kappa))
                s += mp.mpf("0.5")*c[r]*c[t]*(em+ep)
        return s
    cq=corr(coeff,kq);cp=corr(coeff,kp)
    Braw=mp.sqrt(2)*(cq+cp);Bfej=Braw/M8**2
    assert abs(Braw-mp.mpf("2.0949957173756400870"))<mp.mpf("5e-18")
    assert abs(Bfej-mp.mpf("2.2714939727791660463"))<mp.mpf("5e-18")

    Dq=nuA*nuB-kq*kq;dq=nuA+nuB-2*kq
    CR=mp.e**(-mp.pi*dq/(4*Dq))/Dq
    Bhyb=mp.sqrt(2)*(cq/M8**2+CR)
    assert abs(Dq-mp.mpf("1.000493277373620166"))<mp.mpf("5e-18")
    assert abs(dq-mp.mpf("0.02222902596378894775"))<mp.mpf("5e-20")
    assert abs(CR-mp.mpf("0.98221683619439561"))<mp.mpf("5e-17")
    assert abs(Bhyb-mp.mpf("2.55990385767281"))<mp.mpf("5e-14")

    # Optimized hybrid rational coefficients.
    ah=[F("1"),F("-0.312396"),F("0.166775"),F("-0.100361"),
        F("0.061615"),F("-0.036488"),F("0.019464"),F("-0.007871")]
    Mh=F("0.793944")
    P=[F(0)]
    for a,n in zip(ah,ms):P=_padd(P,_pscale(_Tpoly(n),a))
    P2=[F(0)]*(2*len(P)-1)
    for i,a in enumerate(P):
        for j,b in enumerate(P):P2[i+j]+=a*b
    Q=[-x for x in P2];Q[0]+=Mh*Mh
    S=_sturm(Q)
    assert _var([_eval(q,F(-1)) for q in S])==_var([_eval(q,F(1)) for q in S])
    assert _eval(Q,F(-1))>0 and _eval(Q,F(1))>0

    ahmp=[mp.mpf(a.numerator)/a.denominator for a in ah]
    Mact,_=_continuous_norm(ahmp,ms)
    cqhy=corr(ahmp,kq)
    Bcons=mp.sqrt(2)*(cqhy/mp.mpf("0.793944")**2+CR)
    Bact=mp.sqrt(2)*(cqhy/Mact**2+CR)
    assert Bcons>mp.mpf("2.61118871640")
    assert abs(Bact-mp.mpf("2.61118902986463"))<mp.mpf("1e-13")

    print("PASS_SIGMA002_FEJER_REFLECTION_REGRESSIONS")
    print("M8 =",mp.nstr(M8,24))
    print("Fejer Bell =",mp.nstr(Bfej,24))
    print("C_R =",mp.nstr(CR,24))
    print("hybrid Fejer =",mp.nstr(Bhyb,24))
    print("hybrid optimized =",mp.nstr(Bact,24))


def run_sigma001_exact():
    from fractions import Fraction as F
    from math import comb, factorial, isqrt

    # End-to-end exact rational enclosure of the sigma=0.02 Gaussian covariance
    # and degree-15 Bell kernel used by the global optimization certificate.
    SIGMA = F(1,100)
    MS = [1,3,5,7,9,11,13,15]

    class IV:
        def __init__(self, lo, hi=None):
            self.lo=F(lo); self.hi=F(lo if hi is None else hi)
            assert self.lo <= self.hi
        def __add__(self,o):
            o=iv(o); return IV(self.lo+o.lo,self.hi+o.hi)
        __radd__=__add__
        def __neg__(self): return IV(-self.hi,-self.lo)
        def __sub__(self,o): return self+(-iv(o))
        def __rsub__(self,o): return iv(o)-self
        def __mul__(self,o):
            o=iv(o)
            v=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
            return IV(min(v),max(v))
        __rmul__=__mul__
        def recip(self):
            assert not (self.lo <= 0 <= self.hi)
            a,b=F(1,self.lo),F(1,self.hi)
            return IV(min(a,b),max(a,b))
        def __truediv__(self,o): return self*iv(o).recip()
        def __rtruediv__(self,o): return iv(o)/self
        def __pow__(self,n):
            assert isinstance(n,int) and n>=0
            if n==0: return IV(1)
            if self.lo>=0: return IV(self.lo**n,self.hi**n)
            vals=[self.lo**n,self.hi**n]
            if n%2==0 and self.lo<0<self.hi: vals.append(F(0))
            return IV(min(vals),max(vals))
    def iv(x): return x if isinstance(x,IV) else IV(x)

    def atan_inv_bounds(q,N):
        s=F(0)
        for n in range(N):
            t=F(1,(2*n+1)*q**(2*n+1))
            s += t if n%2==0 else -t
        nxt=F(1,(2*N+1)*q**(2*N+1))
        return (s,s+nxt) if N%2==0 else (s-nxt,s)

    # Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    al,ah=atan_inv_bounds(5,80)
    bl,bh=atan_inv_bounds(239,20)
    PI=IV(16*al-4*bh,16*ah-4*bl)

    def sqrt_iv(x,D=80):
        S=10**D
        n=isqrt((x.lo.numerator*S*S)//x.lo.denominator)
        lo=F(n,S)
        while lo*lo>x.lo: n-=1; lo=F(n,S)
        while F(n+1,S)**2<=x.lo: n+=1; lo=F(n,S)
        m=isqrt((x.hi.numerator*S*S)//x.hi.denominator)
        hi=F(m,S)
        if hi*hi<x.hi: m+=1; hi=F(m,S)
        return IV(lo,hi)

    SQRTPI=sqrt_iv(PI,90)

    # Bernoulli numbers and Taylor coefficients.
    B=[F(0)]*31; B[0]=F(1)
    for m in range(1,31):
        B[m]=-sum(F(comb(m+1,k))*B[k] for k in range(m))/F(m+1)
    def cc(m): # z coth z
        return F(1) if m==0 else F(2**(2*m))*B[2*m]/F(factorial(2*m))
    def ct(m): # z tanh z
        return F(0) if m==0 else F(2**(2*m)*(2**(2*m)-1))*B[2*m]/F(factorial(2*m))

    A=PI/200
    PREF=F(1)/(50*PI)

    def moment(m):
        return SQRTPI*F(factorial(2*m),2*(4**m)*factorial(m))

    def partial(kind,M):
        s=IV(0)
        for m in range(M+1):
            c=cc(m) if kind=="A" else ct(m)
            if c: s += c*(A**(2*m))*moment(m)
        return PREF*s

    # On t<=20 the partial-fraction geometric ratios obey q_coth<=0.04,
    # q_tanh<=0.16, so even M=12 is a lower truncation and odd M=13 upper.
    assert F(20*20,10000) <= F(4,100)
    assert F(4*20*20,10000) <= F(16,100)

    # e^-400 < 1e-160 follows from e>2.7 and 2.7^400>1e160.
    assert F(27,10)**400 > 10**160
    E400=F(1,10**160)
    T=F(20)
    def tail_even(m):
        I=E400/(2*T) # I_0
        p=2
        while p<=2*m:
            I=F(1,2)*T**(p-1)*E400+F(p-1,2)*I
            p+=2
        return I
    def poly_tail(kind,M):
        s=F(0)
        for m in range(M+1):
            c=abs(cc(m) if kind=="A" else ct(m))
            if c: s += c*(A.hi**(2*m))*tail_even(m)
        return PREF.hi*s

    I0=tail_even(0); I1=E400/2
    tailA=PREF.hi*(I0+A.hi*I1) # z coth z <= 1+z
    tailB=PREF.hi*(A.hi*I1)    # z tanh z <= z

    A0=IV(partial("A",12).lo-poly_tail("A",12),
          partial("A",13).hi+poly_tail("A",13)+tailA)
    B0=IV(partial("B",12).lo-poly_tail("B",12),
          partial("B",13).hi+poly_tail("B",13)+tailB)

    def exp_pos_small(x,N=14):
        def s(v):
            term=F(1); sm=F(1)
            for n in range(1,N+1):
                term=term*v/F(n); sm+=term
            nxt=term*v/F(N+1)
            return sm,sm+nxt/(1-v/F(N+2))
        lo,_=s(x.lo); _,hi=s(x.hi)
        return IV(lo,hi)

    X=(PI**2)/40000
    EX=exp_pos_small(X,14)
    BC=(PI*SQRTPI)*EX/F(8*10**6)
    AC=A0+BC

    def fourth(x): return sqrt_iv(sqrt_iv(x,75),65)
    NUA=sqrt_iv(AC/BC,65)
    NUB=sqrt_iv(A0/B0,65)
    KQ=fourth(A0**3)/fourth(AC*BC*B0)
    KP=fourth(B0*AC*A0)/fourth(BC**3)

    def round_out(x,N):
        S=10**N
        lo=(x.lo.numerator*S)//x.lo.denominator
        hi=(x.hi.numerator*S+x.hi.denominator-1)//x.hi.denominator
        return IV(F(lo,S),F(hi,S))
    NUA=round_out(NUA,32); NUB=round_out(NUB,32)
    KQ=round_out(KQ,32); KP=round_out(KP,32)

    # Exact e^{-x} bounds for x<=12 by integer+unit-interval alternating series.
    def exp_neg_unit(r,N=70):
        s=F(0); term=F(1)
        for n in range(N+1):
            if n: term*=r/F(n)
            s += term if n%2==0 else -term
        nxt=term*r/F(N+1)
        s2=s+(nxt if (N+1)%2==0 else -nxt)
        return min(s,s2),max(s,s2)
    E1=exp_neg_unit(F(1),70)
    def exp_neg_single(x):
        n=x.numerator//x.denominator; r=x-F(n)
        rl,rh=exp_neg_unit(r,70)
        return E1[0]**n*rl,E1[1]**n*rh
    def exp_neg_iv(x):
        l,_=exp_neg_single(x.hi); _,h=exp_neg_single(x.lo)
        return IV(l,h)

    DQ=NUA+NUB-2*KQ
    DP=NUA+NUB-2*KP

    # The tiny positive-frequency "plus" diagonal terms are bounded separately.
    ELOW=sum(F(1,factorial(n)) for n in range(13)) # strict lower bound for e
    def exp_large_upper(n,r):
        sr=sum(r**k/F(factorial(k)) for k in range(9))
        return F(1)/(ELOW**n*sr)
    PLUS_Q=exp_large_upper(141,F(415,1000))
    PLUS_P=exp_large_upper(141,F(38,100))
    assert (PLUS_Q+PLUS_P)/2 < F(4,10**62)

    KD=[]
    for m in MS:
        xq=round_out((PI/F(4))*F(m*m)*DQ,30)
        xp=round_out((PI/F(4))*F(m*m)*DP,30)
        eq,ep=exp_neg_iv(xq),exp_neg_iv(xp)
        KD.append(IV((eq.lo+ep.lo)/2,(eq.hi+ep.hi)/2+F(4,10**62)))

    # Off-diagonal minus branches: exact lower exponents.
    def qminus(m,n,k): return F(m*m)*NUA+F(n*n)*NUB-F(2*m*n)*k
    minq=min(((PI/F(4))*qminus(m,n,KQ)).lo for i,m in enumerate(MS) for j,n in enumerate(MS) if i!=j)
    minp=min(((PI/F(4))*qminus(m,n,KP)).lo for i,m in enumerate(MS) for j,n in enumerate(MS) if i!=j)
    assert minq > F(141415,1000)
    assert minp > F(14152,100)
    UQ=exp_large_upper(141,F(415,1000))
    UP=exp_large_upper(141,F(52,100))
    # Off-diagonal plus branches have exponent >565, hence <1e-240.
    assert F(27,10)**565 > 10**240
    UOFF=(UQ+UP)/2+F(1,10**240)
    assert UOFF < F(37,10**63) # 3.7e-62
    assert UOFF < F(4,10**62)

    # Exact lower bound for the previously Sturm-certified rational degree-15 witness.
    ahat=[
        F("1.2598224203066"),F("-0.3930260706508"),F("0.2090203393372"),
        F("-0.1250337733856"),F("0.0761080283701"),F("-0.0444368250861"),
        F("0.0231362615839"),F("-0.0089556495026")]
    Mbar=F("1.00000000000006")
    ahat=[a/Mbar for a in ahat]
    qlo=sum(ahat[i]*ahat[i]*KD[i].lo for i in range(8))
    off=F(0)
    for i in range(8):
        for j in range(i+1,8):
            off += 2*abs(ahat[i]*ahat[j])*UOFF
    qlo -= off
    sqrt2_lo=F(1414213562373095,10**15)
    assert sqrt2_lo*sqrt2_lo < 2
    BLOW=sqrt2_lo*qlo
    assert BLOW > F("2.44745254662")

    def dec(x,d=36):
        return f"{float(x):.{d//2}g}" # audit display only

    print("PASS_SIGMA001_WITNESS_END_TO_END_EXACT")
    print("A0 width <", float(A0.hi-A0.lo))
    print("B0 width <", float(B0.hi-B0.lo))
    print("nuA =", float(NUA.lo), float(NUA.hi))
    print("nuB =", float(NUB.lo), float(NUB.hi))
    print("kappa_q =", float(KQ.lo), float(KQ.hi))
    print("kappa_p =", float(KP.lo), float(KP.hi))
    for m,x in zip(MS,KD):
        print("Kdiag",m,float(x.lo),float(x.hi))
    print("offdiag absolute upper <", float(UOFF))
    print("certified witness CHSH lower =", float(BLOW))


def run_bandwidth_continuation():
    """Deterministic numerical continuation for the degree-15 all-plus contact branch.
    The sigma=0.02 and sigma=0.01 rigorous claims are certified separately by exact
    rational replay scripts; this file is the controlled numerical continuation audit.
    """
    import math
    import numpy as np
    from scipy.integrate import quad
    from scipy.optimize import root

    MS=np.arange(1,16,2)

    def coth(x):
        if abs(x)<1e-7:
            return 1/x+x/3-x**3/45
        return 1/math.tanh(x)

    def covariance(sigma):
        def fa(t):
            if t==0: return 2/(math.pi*sigma)
            z=math.pi*sigma*t/2
            return t*coth(z)*math.exp(-t*t)
        def fb(t):
            z=math.pi*sigma*t/2
            return t*math.tanh(z)*math.exp(-t*t)
        A0=sigma*sigma*quad(fa,0,np.inf,epsabs=1e-13,epsrel=1e-13,limit=200)[0]
        B0=sigma*sigma*quad(fb,0,np.inf,epsabs=1e-13,epsrel=1e-13,limit=200)[0]
        Bc=math.pi**1.5*sigma**3/8*math.exp(math.pi**2*sigma**2/4)
        Ac=A0+Bc
        return (math.sqrt(Ac/Bc), math.sqrt(A0/B0),
                A0**0.75/(Ac*Bc*B0)**0.25,
                (B0*Ac*A0)**0.25/Bc**0.75)

    def kernel(sigma):
        if sigma==0: return np.eye(8)
        nuA,nuB,kq,kp=covariance(sigma)
        K=np.zeros((8,8))
        for i,m in enumerate(MS):
            for j,n in enumerate(MS):
                for kap in (kq,kp):
                    xm=math.pi/4*(m*m*nuA+n*n*nuB-2*m*n*kap)
                    xp=math.pi/4*(m*m*nuA+n*n*nuB+2*m*n*kap)
                    K[i,j]+=0.5*math.exp(-xm)+(0 if xp>745 else 0.5*math.exp(-xp))
        return K

    def tv(y):
        T=np.empty(16); T[0]=1.; T[1]=y
        for n in range(1,15): T[n+1]=2*y*T[n]-T[n-1]
        return T[MS]

    def dtv(y):
        out=[]
        for m in MS:
            c=np.zeros(m+1); c[m]=1
            out.append(np.polynomial.chebyshev.chebval(y,np.polynomial.chebyshev.chebder(c)))
        return np.array(out)

    def residual(x,K):
        a=x[:8]; lam=x[8:12]; y=x[12:]
        V=np.array([tv(z) for z in y]); D=np.array([dtv(z) for z in y])
        return np.r_[V@a-np.ones(4),D@a,K@a-V.T@lam]

    x=np.r_[
     [1.2598224203066,-0.3930260706508,0.2090203393372,-0.1250337733856,
      0.0761080283701,-0.0444368250861,0.0231362615839,-0.0089556495026],
     [0.31798191443,0.43452134604,0.45769155766,0.46184185531],
     [0.319388191001,0.620152408461,0.855496118401,0.983520510350]]

    schedule=[.02,.015,.01,.005,.002,.001,.0005,.0001,0.0]
    for s in schedule:
        K=kernel(s)
        sol=root(lambda z: residual(z,K),x,tol=1e-12)
        if not sol.success or np.linalg.norm(residual(sol.x,K))>1e-8:
            raise RuntimeError((s,sol.message,np.linalg.norm(residual(sol.x,K))))
        x=sol.x
        B=math.sqrt(2)*x[:8]@K@x[:8]
        print(f"{s:.7g}  B={B:.15f}  contacts={x[12:]}")

    # First-order zero-bandwidth envelope coefficient on the continued branch.
    a=x[:8]
    s2=float(np.sum((MS**2)*(a**2)))
    c1=math.pi**2*s2/4
    print("zero-bandwidth candidate =",math.sqrt(2)*float(a@a))
    print("first-order coefficient =",c1)
    assert abs(B-2.56173695411067) < 5e-11
    assert abs(c1-14.81233586) < 5e-6
    print("PASS_BANDWIDTH_CONTINUATION_REGRESSION")

def main():
    run_sigma002_regressions()
    run_sigma001_exact()
    run_bandwidth_continuation()
    print("PASS_GITHUB_TEST_03_BANDWIDTH_REFLECTION_REGRESSIONS")

def test_all_bandwidth_reflection():
    main()

if __name__=="__main__":
    main()
