from sympy import *
from ga import *

#coords = (x, y) = symbols('x y')
(ex, ey) = MV.setup('e*x|y',metric='[1,1]')


# (a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) = symbols('a b c d e f g h i j k l m n o p')
# (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P) = symbols('A B C D E F G H I J K L M N O P')

# MV1 = GA2E(a, b,c, d);
(a1, b1, c1, d1) = symbols('a.q   a.x a.y  a.xy')

(a2, b2, c2, d2) = symbols('b.q   b.x b.y  b.xy')

exy = ex*ey

# ******************************************************

# Test the generic MV product using default basis

print
print '**************************************************'
print
print 'Metric + + Geometric Product'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = a2 + b2*ex + c2*ey + d2*exy
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
print 'Metric + + Inner Product'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = a2 + b2*ex + c2*ey + d2*exy
BMV.Fmt(2,'BMV')

CMV = (AMV|BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Metric + + Left Contraction Product'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = a2 + b2*ex + c2*ey + d2*exy
BMV.Fmt(2,'BMV')


CMV = (AMV<BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Metric + + Right Contraction Product'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = a2 + b2*ex + c2*ey + d2*exy
BMV.Fmt(2,'BMV')

CMV = (AMV>BMV)
CMV.Fmt(3,'CMV')
print


print
print '**************************************************'
print
print 'Dual = AMV*(I_inv) = AMV*exy'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = AMV*(-exy)
BMV.Fmt(2,'BMV')

print


print
print '**************************************************'
print
print '	DorstDual b = LeftContraction(a,I_inv);'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = AMV<(-exy)
BMV.Fmt(2,'BMV')

print


print
print '**************************************************'
print
print '	DorstUnDual b = LeftContraction(a,I);'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = AMV<(exy)
BMV.Fmt(2,'BMV')

print


print
print '**************************************************'
print
print '	Hestenes Regressive = c = ((a*I_inv)^(b*I_inv))*I;'
print 

AMV = a1 + b1*ex + c1*ey + d1*exy
AMV.Fmt(2,'AMV')

BMV = a2 + b2*ex + c2*ey + d2*exy
BMV.Fmt(2,'BMV')

CMV = ((AMV*(-exy))^(BMV*(-exy)))*exy
CMV.Fmt(3,'CMV')
print








