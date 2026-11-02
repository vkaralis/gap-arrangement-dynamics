"""Exact rational bound for the affine observation-time example, not T=0.02."""
from fractions import Fraction as F


def certificate():
    r=F(22,7)
    K0,K1,K2=F(3,2),r,2*r*r
    w1,w2,w3=2*r,(2*r)**2,(2*r)**3
    B,B2=12*r*r,32*r**3
    S0=F(15,7)
    S1,S2,S3,S4=[F(3,7)*(4*r)**j for j in range(1,5)]
    M0,M1=1+K0,K1
    V0=M0*w1
    Q1=B*S0+S2
    Q11=B*S1+S3
    Q12=B*S2+S4
    R0=K1*V0+K0*Q1
    R1=K2*V0+K1*Q1
    T0=R0*w1+M0*w2*V0
    T1=R1*w1+R0*w2+M1*w2*V0+M0*w3*V0
    Q2=(T1+B2*V0)*S0+B*Q1+Q12
    V2=T0+B*V0
    Rtt=K2*V0**2+K1*V2+2*K1*V0*Q1+K0*(Q1**2+Q2)
    T2=Rtt*w1+2*R0*w2*V0+M0*(w3*V0**2+w2*V2)
    V3=T2+2*T1*V0+B2*V0**2+B*V2
    L3=2*(w3*V0**3+3*w2*V0*V2+w1*V3)
    cap=10**len(str(L3.numerator//L3.denominator))
    assert L3<cap
    tcert=3*F(243,4096)/(2*cap)
    assert cap == 10**9
    assert tcert == F(729,8192*10**9)
    assert 12*r*r/F(1000) < F(3,25)
    assert F(25,22) < F(8,7)
    return {
        'L3_upper_evaluation': float(L3),
        'L3_upper_rational': str(L3),
        'cap': cap,
        'certified_time': float(tcert),
        'certified_time_rational': str(tcert),
        'scope': 'Affine theoretical example only; not the simulation horizon T=0.02.',
    }
