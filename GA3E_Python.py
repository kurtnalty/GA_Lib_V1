from sympy import *
from ga import *

#coords = (x, y, z) = symbols('x y z')
(ex, ey, ez) = MV.setup('e*x|y|z',metric='[1,1,1]')


# (a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) = symbols('a b c d e f g h i j k l m n o p')
# (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P) = symbols('A B C D E F G H I J K L M N O P')

(a, b, c, d, e, f, g, h) = symbols('a.q   a.x a.y a.z   a.xy a.xz a.yz    a.xyz')
(A, B, C, D, E, F, G, H) = symbols('b.q   b.x b.y b.z   b.xy b.xz b.yz    b.xyz')


# ******************************************************

# Print the generic MV product using default basis

print
print '**************************************************'
print
print 'Euclidean 3D Metric + + + Geometric Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*ex*ey + F*ex*ez + G*ey*ez + H*ex*ey*ez
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
print 'Euclidean 3D Metric + + + Inner Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*ex*ey + F*ex*ez + G*ey*ez + H*ex*ey*ez
BMV.Fmt(2,'BMV')

CMV = (AMV|BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Euclidean 3D Metric + + + Left Contraction Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*ex*ey + F*ex*ez + G*ey*ez + H*ex*ey*ez
BMV.Fmt(2,'BMV')

CMV = (AMV<BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Euclidean 3D Metric + + + Right Contraction Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*ex*ey + F*ex*ez + G*ey*ez + H*ex*ey*ez
BMV.Fmt(2,'BMV')

CMV = (AMV>BMV)
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Euclidean 3D Metric + + + Symmetric Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*ex*ey + F*ex*ez + G*ey*ez + H*ex*ey*ez
BMV.Fmt(2,'BMV')

CMV = ((AMV*BMV) + (BMV*AMV))/2
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Euclidean 3D Metric + + + AntiSymmetric Product'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*ex*ey + F*ex*ez + G*ey*ez + H*ex*ey*ez
BMV.Fmt(2,'BMV')

CMV = ((AMV*BMV) - (BMV*AMV))/2
CMV.Fmt(3,'CMV')
print

print
print '**************************************************'
print
print 'Euclidean 3D Metric + + + Dual = w*I_inv'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = - ex*ey*ez
BMV.Fmt(2,'BMV')

CMV = AMV*BMV
CMV.Fmt(3,'CMV')
print


print
print '**************************************************'
print
print 'Euclidean 3D Metric + + + Regressive = ((a*I_inv)^(b*I_inv))*I'
print 
AMV = a + b*ex + c*ey + d*ez + e*ex*ey + f*ex*ez + g*ey*ez + h*ex*ey*ez
AMV.Fmt(2,'AMV')

BMV = A + B*ex + C*ey + D*ez + E*ex*ey + F*ex*ez + G*ey*ez + H*ex*ey*ez
BMV.Fmt(2,'BMV')

CMV = ((AMV*(-ex*ey*ez))^(BMV*(-ex*ey*ez)) )*ex*ey*ez
CMV.Fmt(3,'CMV')
print







