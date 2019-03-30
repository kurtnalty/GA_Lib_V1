from sympy import *
from ga import *

#coords = (x, y, z, t) = symbols('x y z t')
(ex, ey, ez, et) = MV.setup('e*x|y|z|t',metric='[1,1,1,1]')


# (a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) = symbols('a b c d e f g h i j k l m n o p')
# (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P) = symbols('A B C D E F G H I J K L M N O P')

(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) = symbols('a.q   a.x a.y a.z a.t   a.xy a.xz a.xt a.yz a.yt a.zt    a.xyz a.xyt a.xzt a.yzt    a.xyzt')
(A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P) = symbols('b.q   b.x b.y b.z b.t   b.xy b.xz b.xt b.yz b.yt b.zt    b.xyz b.xyt b.xzt b.yzt    b.xyzt')


# ******************************************************

# Test the generic MV product using default basis

print
print '**************************************************'
print
print 'Euclidean 4D Metric + + + + Geometric Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*et + f*ex*ey + g*ex*ez + h*ex*et + i*ey*ez + j*ey*et + k*ez*et + l*ex*ey*ez + m*ex*ey*et + n*ex*ez*et + o*ey*ez*et + p*ex*ey*ez*et
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*et + F*ex*ey + G*ex*ez + H*ex*et + I*ey*ez + J*ey*et + K*ez*et + L*ex*ey*ez + M*ex*ey*et + N*ex*ez*et + O*ey*ez*et + P*ex*ey*ez*et
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
print 'Euclidean 4D Metric + + + + Inner Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*et + f*ex*ey + g*ex*ez + h*ex*et + i*ey*ez + j*ey*et + k*ez*et + l*ex*ey*ez + m*ex*ey*et + n*ex*ez*et + o*ey*ez*et + p*ex*ey*ez*et
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*et + F*ex*ey + G*ex*ez + H*ex*et + I*ey*ez + J*ey*et + K*ez*et + L*ex*ey*ez + M*ex*ey*et + N*ex*ez*et + O*ey*ez*et + P*ex*ey*ez*et
BMV.Fmt(2,'BMV')

CMV = (AMV|BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Euclidean 4D Metric + + + + Left Contraction Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*et + f*ex*ey + g*ex*ez + h*ex*et + i*ey*ez + j*ey*et + k*ez*et + l*ex*ey*ez + m*ex*ey*et + n*ex*ez*et + o*ey*ez*et + p*ex*ey*ez*et
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*et + F*ex*ey + G*ex*ez + H*ex*et + I*ey*ez + J*ey*et + K*ez*et + L*ex*ey*ez + M*ex*ey*et + N*ex*ez*et + O*ey*ez*et + P*ex*ey*ez*et
BMV.Fmt(2,'BMV')

CMV = (AMV<BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Euclidean 4D Metric + + + + Right Contraction Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*et + f*ex*ey + g*ex*ez + h*ex*et + i*ey*ez + j*ey*et + k*ez*et + l*ex*ey*ez + m*ex*ey*et + n*ex*ez*et + o*ey*ez*et + p*ex*ey*ez*et
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*et + F*ex*ey + G*ex*ez + H*ex*et + I*ey*ez + J*ey*et + K*ez*et + L*ex*ey*ez + M*ex*ey*et + N*ex*ez*et + O*ey*ez*et + P*ex*ey*ez*et
BMV.Fmt(2,'BMV')

CMV = (AMV>BMV)
CMV.Fmt(3,'CMV')
print








