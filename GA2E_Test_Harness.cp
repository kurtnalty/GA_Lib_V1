#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

#include "GA2E_Routines.cp"

/////////////////////////////////// Demonstrate GA2E Routines ////////////////////////////////


int main(void)
{
	
	GA2E Zero;   Zero = GA2E();

	symbol a_q("a"), a_x("b"), a_y("c"), a_xy("d") ;
	symbol b_q("A"), b_x("B"), b_y("C"), b_xy("D") ;
	symbol c_q("c.a"), c_x("c.b"), c_y("c.c"), c_xy("c.d") ;
	symbol d_q("d.a"), d_x("d.b"), d_y("d.c"), d_xy("d.d") ;

	GA2E r,s,t,u,w;

	GA2E X, Y, Z;

	X = GA2E( 3,  5,  7, 11) ; 
	Y = GA2E(13, 17, 19, 23) ;
	Z = GA2E(29, 31, 37, 41) ; 

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

	ex det_ref, det;
	det_ref = Determinant(r);

//////////////////////////////////////////////////////

// ostream &operator<<(ostream &ff, GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	cout << "\nr = " << r << "\n";

//////////////////////////////////////////////////////

// void PrintMV(GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	cout << "r = "; PrintMV(r); cout << "\n";

//////////////////////////////////////////////////////

// GA2E OverBar(GA2E a) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = OverBar(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "OverBar(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*OverBar(r) = " << s  << "   expect zero for s.q\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E UnderBar(GA2E a) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	u = UnderBar(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "UnderBar(r) = "; PrintMV(u); cout << "\n";
	s = OverBar(u);
	cout << "OverBar(UnderBar(r)) = "; PrintMV(s); cout << "\n";
	u = OverBar(r);	s = UnderBar(u);
	cout << "UnderBar(OverBar(r)) = "; PrintMV(s); cout << "\n";

	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";

	r = GA2E(a_q, a_x,a_y, a_xy);  
	u = UnderBar(r);
	s = r*u;
	cout << "r*UnderBar(r) = " << s  << "   expect zeroes in s.q, s.x, s.y\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E Reverse(GA2E w) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	u = Reverse(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "Reverse(r) = "; PrintMV(u); cout << "\n";
	s = Reverse(u);
	cout << "Reverse(Reverse(r)) = "; PrintMV(s); cout << "\n";

	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*Reverse(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E Involution(GA2E w) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = Involution(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "Involution(r) = "; PrintMV(u); cout << "\n";
	s = Involution(u);
	cout << "Involution(Involution(r)) = "; PrintMV(s); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*Involution(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E Transpose(GA2E w) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = Transpose(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "Transpose(r) = "; PrintMV(u); cout << "\n";
	s = Transpose(u);
	cout << "Transpose(Transpose(r)) = "; PrintMV(s); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*Transpose(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E Conjugation(GA2E w) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = Conjugation(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "Conjugation(r) = "; PrintMV(u); cout << "\n";
	s = Conjugation(u);
	cout << "Conjugation(Conjugation(r)) = "; PrintMV(s); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*Conjugation(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E CliffordConjugation(GA2E w) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = CliffordConjugation(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "CliffordConjugation(r) = "; PrintMV(u); cout << "\n";
	s = CliffordConjugation(u);
	cout << "CliffordConjugation(CliffordConjugation(r)) = "; PrintMV(s); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*CliffordConjugation(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E Dual(GA2E w) ; 

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = Dual(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "Dual(r) = "; PrintMV(u); cout << "\n";
	s = Dual(u);
	cout << "Dual(Dual(r)) = "; PrintMV(s); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*Dual(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E DorstDual(GA2E a) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = DorstDual(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "DorstDual(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*DorstDual(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E DorstUnDual(GA2E a) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  

	u = DorstUnDual(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "DorstUnDual(r) = "; PrintMV(u); cout << "\n";
	cout << "\n";

	s = DorstDual(r);
	u = DorstUnDual(s);

	cout << "DorstUnDual(DorstDual(r)) = "; PrintMV(u); cout << "\n";

	s = DorstUnDual(r);
	u = DorstDual(s);

	cout << "DorstDual(DorstUnDual(r)) = "; PrintMV(u); cout << "\n";

	u = DorstUnDual(r);
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*DorstUnDual(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E operator+(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 

	u = r + s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r + s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E operator-(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	u = r - s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r - s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator==(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(b_q, b_x,b_y, b_xy); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r == s) cout << "r == s\n";   else cout << "r != s (as expected)\n";
	if (s == t) cout << "s == t (as expected)\n";   else cout << "s != t \n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator!=(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(b_q, b_x,b_y, b_xy); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r != s) cout << "r != s (as expected)\n";   else cout << "r == s \n";
	if (s != t) cout << "s != t \n";   else cout << "s == t (as expected)\n";
	cout << "\n";

//////////////////////////////////////////////////////

GA2E operator*(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r*s;
	cout << "u = r*s = " << u << "\n";

	u = r*s - s*r;
	if (u == Zero) cout << "Product is commutative\n"; else cout << "Product is non-commutative (as expected)\n";
	cout << "\n";

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = (r*s)*t - r*(s*t);
	if (u == Zero) cout << "Product is associative (as expected)\n"; else cout << "Product is non-associative\n";
	cout << "\n";

	s = symbol("Z")*r;
	cout << "s = Z*r = "; PrintMV(s); cout << "\n";
	s = r*symbol("Z");
	cout << "s = r*Z = "; PrintMV(s); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E Product(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Product(r,s);
	cout << "u = Product(r,s) = " << u << "\n";

	u = Product(r,s) - Product(s,r);
	if (u == Zero) cout << "Product(r,s) is commutative\n"; else cout << "Product(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Product(Product(r,s), t) - Product(r,Product(s,t));
	if (u == Zero) cout << "Product is associative (as expected)\n"; else cout << "Product is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E operator/(const GA2E &u, const int i) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	cout << "r = "; PrintMV(r); cout << "\n";
	u = r/2;
	cout << "u = r/2 = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E operator^(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r^s;
	cout << "u = r^s = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E Wedge(const GA2E &u, const GA2E &v) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Wedge(r,s);
	cout << "u = Wedge(r,s) = " << u << "\n";

	u = Wedge(r,s) - Wedge(s,r);
	if (u == Zero) cout << "Wedge(r,s) is commutative\n"; else cout << "Wedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Wedge(Wedge(r,s), t) - Wedge(r,Wedge(s,t));
	if (u == Zero) cout << "Wedge is associative (as expected)\n"; else cout << "Wedge is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E AntiWedge(const GA2E a, const GA2E b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiWedge(r,s);
	cout << "u = AntiWedge(r,s) = " << u << "\n";

	u = AntiWedge(r,s) - AntiWedge(s,r);
	if (u == Zero) cout << "AntiWedge(r,s) is commutative\n"; else cout << "AntiWedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = AntiWedge(AntiWedge(r,s), t) - AntiWedge(r,AntiWedge(s,t));
	if (u == Zero) cout << "AntiWedge is associative (as expected)\n"; else cout << "AntiWedge is non-associative\n";
	cout << "\n";


//////////////////////////////////////////////////////

// GA2E Regressive(GA2E a, GA2E b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Regressive(r,s);
	cout << "u = Regressive(r,s) = " << u << "\n";

	u = Regressive(r,s) - Regressive(s,r);
	if (u == Zero) cout << "Regressive(r,s) is commutative\n"; else cout << "Regressive(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Regressive(Regressive(r,s), t) - Regressive(r,Regressive(s,t));
	if (u == Zero) cout << "Regressive is associative (as expected)\n"; else cout << "Regressive is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA2E RegressiveViaFormula(GA2E a, GA2E b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RegressiveViaFormula(r,s);
	cout << "u = RegressiveViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E LowerRightViaFormula(GA2E a, GA2E b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LowerRightViaFormula(r,s);
	cout << "u = LowerRightViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E Expander(const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Expander(r,s);
	cout << "u = Expander(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E Conserver(const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Conserver(r,s);
	cout << "u = Conserver(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E Shrinker(const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Shrinker(r,s);
	cout << "u = Shrinker(r,s) = " << u << "\n\n";

///////////////////////////////////////////////////////

// GA2E Symmetric(const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Symmetric(r,s);
	cout << "u = Symmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E AntiSymmetric(const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiSymmetric(r,s);
	cout << "u = AntiSymmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E Inner(const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Inner(r,s);
	cout << "u = Inner(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E LeftContraction (const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LeftContraction(r,s);
	cout << "u = LeftContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E RightContraction (const GA2E &a, const GA2E &b) ;

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RightContraction(r,s);
	cout << "u = RightContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// ex Determinant(GA2E A) ;

	X = GA2E( 3,  5,  7, 11) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	cout << "Determinant(X) = " << Determinant(X) << "\n\n";


//////////////////////////////////////////////////////

// GA2E Adjugate(GA2E V) ;

	X = GA2E( 3,  5,  7, 11) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Adjugate(X);
	cout << "u = Adjugate(X) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA2E Reciprocal(GA2E a) ;

	X = GA2E( 3,  5,  7, 11) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Reciprocal(X);
	cout << "u = Reciprocal(X) = " << u << "\n";
	s = X*u;
	cout << "X*(1/X) = " << s << "\n\n";

/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	X = GA2E( 3,  5,  7, 11) ; 
	Y = GA2E(13, 17, 19, 23) ;
	Z = GA2E(29, 31, 37, 41) ; 

	r = GA2E(a_q, a_x,a_y, a_xy);  
	s = GA2E(b_q, b_x,b_y, b_xy); 
	t = GA2E(c_q, c_x,c_y, c_xy); 
	u = GA2E(d_q, d_x,d_y, d_xy); 

*/





