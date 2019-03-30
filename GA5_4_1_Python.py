from sympy import *
from ga import *

#coords = (w, x, y, z, t) = symbols('w x y z t')
(ew, ex, ey, ez, et) = MV.setup('e*w|x|y|z|t',metric='[1,1,1,1,-1]')


# (a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) = symbols('a b c d e f g h i j k l m n o p')
# (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P) = symbols('A B C D E F G H I J K L M N O P')

# MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
(a1, b1, c1, d1, e1, f1) = symbols('a.q   a.w a.x a.y a.z a.t')
(g1, h1, j1, k1, l1, m1, n1, p1, r1, s1) = symbols('a.wx  a.wy  a.wz  a.wt  a.xy  a.xz  a.xt  a.yz  a.yt  a.zt ')   
(S1, R1, P1, N1, M1, L1, K1, J1, H1, G1) = symbols('a.wxy a.wxz a.wxt a.wyz a.wyt a.wzt a.xyz a.xyt a.xzt a.yzt')
(F1, E1, D1, C1, B1, A1) = symbols('a.wxyz a.wxyt a.wxzt a.wyzt a.xyzt a.wxyzt')

(a2, b2, c2, d2, e2, f2) = symbols('b.q   b.w b.x b.y b.z b.t')
(g2, h2, j2, k2, l2, m2, n2, p2, r2, s2) = symbols('b.wx  b.wy  b.wz  b.wt  b.xy  b.xz  b.xt  b.yz  b.yt  b.zt ')   
(S2, R2, P2, N2, M2, L2, K2, J2, H2, G2) = symbols('b.wxy b.wxz b.wxt b.wyz b.wyt b.wzt b.xyz b.xyt b.xzt b.yzt')
(F2, E2, D2, C2, B2, A2) = symbols('b.wxyz b.wxyt b.wxzt b.wyzt b.xyzt b.wxyzt')

ewx = ew*ex
ewy = ew*ey
ewz = ew*ez
ewt = ew*et
exy = ex*ey
exz = ex*ez
ext = ex*et
eyz = ey*ez
eyt = ey*et
ezt = ez*et

ewxy = ewx*ey
ewxz = ewx*ez
ewxt = ewx*et
ewyz = ewy*ez
ewyt = ewy*et
ewzt = ewz*et
exyz = exy*ez
exyt = exy*et
exzt = exz*et
eyzt = eyz*et

ewxyz = ewxy*ez
ewxyt = ewxy*et
ewxzt = ewxz*et
ewyzt = ewyz*et
exyzt = exyz*et
ewxyzt = ewxyz*et

# ******************************************************

# Test the generic MV product using default basis

print
print '**************************************************'
print
print 'Metric + + + + - Geometric Product'
print 

AM1 = a1 + b1*ew + c1*ex + d1*ey + e1*ez + f1*et
AM2 = g1*ewx  + h1*ewy +  j1*ewz +  k1*ewt +  l1*exy +  m1*exz +  n1*ext +  p1*eyz +  r1*eyt +  s1*ezt
AM3 = S1*ewxy + R1*ewxz + P1*ewxt + N1*ewyz + M1*ewyt + L1*ewzt + K1*exyz + J1*exyt + H1*exzt + G1*eyzt
AM4 = F1*ewxyz + E1*ewxyt + D1*ewxzt + C1*ewyzt + B1*exyzt + A1*ewxyzt
AMV = AM1 + AM2 + AM3 + AM4

AMV.Fmt(2,'AMV')

BM1 = a2 + b2*ew + c2*ex + d2*ey + e2*ez + f2*et
BM2 = g2*ewx  + h2*ewy +  j2*ewz +  k2*ewt +  l2*exy +  m2*exz +  n2*ext +  p2*eyz +  r2*eyt +  s2*ezt
BM3 = S2*ewxy + R2*ewxz + P2*ewxt + N2*ewyz + M2*ewyt + L2*ewzt + K2*exyz + J2*exyt + H2*exzt + G2*eyzt
BM4 = F2*ewxyz + E2*ewxyt + D2*ewxzt + C2*ewyzt + B2*exyzt + A2*ewxyzt
BMV = BM1 + BM2 + BM3 + BM4

BMV.Fmt(2,'BMV')

CMV = AMV*BMV
CMV.Fmt(3,'CMV')
print

# ******************************************************

print
print '**************************************************'
print
print 'Wedge Product'
print 

CMV = (AMV^BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Metric + + + + - Inner Product'
print 

AM1 = a1 + b1*ew + c1*ex + d1*ey + e1*ez + f1*et
AM2 = g1*ewx  + h1*ewy +  j1*ewz +  k1*ewt +  l1*exy +  m1*exz +  n1*ext +  p1*eyz +  r1*eyt +  s1*ezt
AM3 = S1*ewxy + R1*ewxz + P1*ewxt + N1*ewyz + M1*ewyt + L1*ewzt + K1*exyz + J1*exyt + H1*exzt + G1*eyzt
AM4 = F1*ewxyz + E1*ewxyt + D1*ewxzt + C1*ewyzt + B1*exyzt + A1*ewxyzt
AMV = AM1 + AM2 + AM3 + AM4

AMV.Fmt(2,'AMV')

BM1 = a2 + b2*ew + c2*ex + d2*ey + e2*ez + f2*et
BM2 = g2*ewx  + h2*ewy +  j2*ewz +  k2*ewt +  l2*exy +  m2*exz +  n2*ext +  p2*eyz +  r2*eyt +  s2*ezt
BM3 = S2*ewxy + R2*ewxz + P2*ewxt + N2*ewyz + M2*ewyt + L2*ewzt + K2*exyz + J2*exyt + H2*exzt + G2*eyzt
BM4 = F2*ewxyz + E2*ewxyt + D2*ewxzt + C2*ewyzt + B2*exyzt + A2*ewxyzt
BMV = BM1 + BM2 + BM3 + BM4

BMV.Fmt(2,'BMV')

CMV = (AMV|BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Metric + + + + - Left Contraction Product'
print 

AM1 = a1 + b1*ew + c1*ex + d1*ey + e1*ez + f1*et
AM2 = g1*ewx  + h1*ewy +  j1*ewz +  k1*ewt +  l1*exy +  m1*exz +  n1*ext +  p1*eyz +  r1*eyt +  s1*ezt
AM3 = S1*ewxy + R1*ewxz + P1*ewxt + N1*ewyz + M1*ewyt + L1*ewzt + K1*exyz + J1*exyt + H1*exzt + G1*eyzt
AM4 = F1*ewxyz + E1*ewxyt + D1*ewxzt + C1*ewyzt + B1*exyzt + A1*ewxyzt
AMV = AM1 + AM2 + AM3 + AM4

AMV.Fmt(2,'AMV')

BM1 = a2 + b2*ew + c2*ex + d2*ey + e2*ez + f2*et
BM2 = g2*ewx  + h2*ewy +  j2*ewz +  k2*ewt +  l2*exy +  m2*exz +  n2*ext +  p2*eyz +  r2*eyt +  s2*ezt
BM3 = S2*ewxy + R2*ewxz + P2*ewxt + N2*ewyz + M2*ewyt + L2*ewzt + K2*exyz + J2*exyt + H2*exzt + G2*eyzt
BM4 = F2*ewxyz + E2*ewxyt + D2*ewxzt + C2*ewyzt + B2*exyzt + A2*ewxyzt
BMV = BM1 + BM2 + BM3 + BM4

BMV.Fmt(2,'BMV')

CMV = (AMV<BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Metric + + + + - Right Contraction Product'
print 
AM1 = a1 + b1*ew + c1*ex + d1*ey + e1*ez + f1*et
AM2 = g1*ewx  + h1*ewy +  j1*ewz +  k1*ewt +  l1*exy +  m1*exz +  n1*ext +  p1*eyz +  r1*eyt +  s1*ezt
AM3 = S1*ewxy + R1*ewxz + P1*ewxt + N1*ewyz + M1*ewyt + L1*ewzt + K1*exyz + J1*exyt + H1*exzt + G1*eyzt
AM4 = F1*ewxyz + E1*ewxyt + D1*ewxzt + C1*ewyzt + B1*exyzt + A1*ewxyzt
AMV = AM1 + AM2 + AM3 + AM4

AMV.Fmt(2,'AMV')

BM1 = a2 + b2*ew + c2*ex + d2*ey + e2*ez + f2*et
BM2 = g2*ewx  + h2*ewy +  j2*ewz +  k2*ewt +  l2*exy +  m2*exz +  n2*ext +  p2*eyz +  r2*eyt +  s2*ezt
BM3 = S2*ewxy + R2*ewxz + P2*ewxt + N2*ewyz + M2*ewyt + L2*ewzt + K2*exyz + J2*exyt + H2*exzt + G2*eyzt
BM4 = F2*ewxyz + E2*ewxyt + D2*ewxzt + C2*ewyzt + B2*exyzt + A2*ewxyzt
BMV = BM1 + BM2 + BM3 + BM4

BMV.Fmt(2,'BMV')

CMV = (AMV>BMV)
CMV.Fmt(3,'CMV')
print








