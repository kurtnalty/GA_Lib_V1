#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

#include "GA3E_Routines.cp"

/////////////////////////////////// Demonstrate GA3E Routines ////////////////////////////////


int main(void)
{
	
	GA3E Zero;   Zero = GA3E();

	symbol a_q("a");
	symbol a_x("b") ,  a_y("c") ,  a_z("d");
	symbol a_xy("e"),  a_xz("f"),  a_yz("g");
	symbol a_xyz("h");

	symbol b_q("A");
	symbol b_x("B") ,  b_y("C") ,  b_z("D");
	symbol b_xy("E"),  b_xz("F"),  b_yz("G");
	symbol b_xyz("H");

	symbol c_q("c.a");
	symbol c_x("c.b") ,  c_y("c.c") ,  c_z("c.d"),   c_t("c.e");
	symbol c_xy("c.f"),  c_xz("c.g"),  c_yz("c.h"),  c_xt("c.i"), c_yt("c.j"), c_zt("c.k");
	symbol c_xyz("c.l"), c_xyt("c.m"), c_xzt("c.n"), c_yzt("c.o");
	symbol c_xyzt("c.p");

	symbol d_q("d.a");
	symbol d_x("d.b") ,  d_y("d.c") ,  d_z("d.d"),   d_t("d.e");
	symbol d_xy("d.f"),  d_xz("d.g"),  d_yz("d.h"),  d_xt("d.i"), d_yt("d.j"), d_zt("d.k");
	symbol d_xyz("d.l"), d_xyt("d.m"), d_xzt("d.n"), d_yzt("d.o");
	symbol d_xyzt("d.p");

	GA3E r,s,t,u,w;

	GA3E X, Y, Z;

	X = GA3E(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	Y = GA3E( 29,   31, 37, 41,     43, 47, 53,   59) ;
	Z = GA3E( 61,   67, 71, 73,     79, 83, 89,   97) ; 

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	u = GA3E(d_q, d_x,d_y,d_z, d_xy,d_xz,d_yz, d_xyz); 

	ex det_ref, det;
	det_ref = Determinant(r);

//////////////////////////////////////////////////////

// ostream &operator<<(ostream &ff, GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	cout << "\nr = " << r << "\n";

//////////////////////////////////////////////////////

// void PrintMV(GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	cout << "r = "; PrintMV(r); cout << "\n";

//////////////////////////////////////////////////////

// GA3E OverBar(GA3E a) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

	u = OverBar(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "OverBar(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*OverBar(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E UnderBar(GA3E a) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

	u = UnderBar(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "UnderBar(r) = "; PrintMV(u); cout << "\n";
	s = OverBar(u);
	cout << "OverBar(UnderBar(r)) = "; PrintMV(s); cout << "\n";
	u = OverBar(r);	s = UnderBar(u);
	cout << "UnderBar(OverBar(r)) = "; PrintMV(s); cout << "\n";

	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*UnderBar(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E Reverse(GA3E w) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

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

// GA3E Involution(GA3E w) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

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

// GA3E Transpose(GA3E w) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

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

// GA3E Conjugation(GA3E w) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

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

// GA3E CliffordConjugation(GA3E w) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

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

// GA3E Dual(GA3E w) ; 

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

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

// GA3E DorstDual(GA3E a) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

	u = DorstDual(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "DorstDual(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*DorstDual(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E DorstUnDual(GA3E a) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  

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

// GA3E operator+(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	u = r + s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r + s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E operator-(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	u = r - s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r - s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator==(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r == s) cout << "r == s\n";   else cout << "r != s (as expected)\n";
	if (s == t) cout << "s == t (as expected)\n";   else cout << "s != t \n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator!=(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r != s) cout << "r != s (as expected)\n";   else cout << "r == s \n";
	if (s != t) cout << "s != t \n";   else cout << "s == t (as expected)\n";
	cout << "\n";

//////////////////////////////////////////////////////

GA3E operator*(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r*s;
	cout << "u = r*s = " << u << "\n";

	u = r*s - s*r;
	if (u == Zero) cout << "Product is commutative\n"; else cout << "Product is non-commutative (as expected)\n";
	cout << "\n";

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
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

// GA3E Product(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Product(r,s);
	cout << "u = Product(r,s) = " << u << "\n";

	u = Product(r,s) - Product(s,r);
	if (u == Zero) cout << "Product(r,s) is commutative\n"; else cout << "Product(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Product(Product(r,s), t) - Product(r,Product(s,t));
	if (u == Zero) cout << "Product is associative (as expected)\n"; else cout << "Product is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E operator/(const GA3E &u, const int i) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	cout << "r = "; PrintMV(r); cout << "\n";
	u = r/2;
	cout << "u = r/2 = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E operator^(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r^s;
	cout << "u = r^s = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E Wedge(const GA3E &u, const GA3E &v) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Wedge(r,s);
	cout << "u = Wedge(r,s) = " << u << "\n";

	u = Wedge(r,s) - Wedge(s,r);
	if (u == Zero) cout << "Wedge(r,s) is commutative\n"; else cout << "Wedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Wedge(Wedge(r,s), t) - Wedge(r,Wedge(s,t));
	if (u == Zero) cout << "Wedge is associative (as expected)\n"; else cout << "Wedge is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E AntiWedge(const GA3E a, const GA3E b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiWedge(r,s);
	cout << "u = AntiWedge(r,s) = " << u << "\n";

	u = AntiWedge(r,s) - AntiWedge(s,r);
	if (u == Zero) cout << "AntiWedge(r,s) is commutative\n"; else cout << "AntiWedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = AntiWedge(AntiWedge(r,s), t) - AntiWedge(r,AntiWedge(s,t));
	if (u == Zero) cout << "AntiWedge is associative (as expected)\n"; else cout << "AntiWedge is non-associative\n";
	cout << "\n";


//////////////////////////////////////////////////////

// GA3E Regressive(GA3E a, GA3E b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Regressive(r,s);
	cout << "u = Regressive(r,s) = " << u << "\n";

	u = Regressive(r,s) - Regressive(s,r);
	if (u == Zero) cout << "Regressive(r,s) is commutative\n"; else cout << "Regressive(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	t = GA3E(c_q, c_x,c_y,c_z, c_xy,c_xz,c_yz, c_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Regressive(Regressive(r,s), t) - Regressive(r,Regressive(s,t));
	if (u == Zero) cout << "Regressive is associative (as expected)\n"; else cout << "Regressive is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA3E RegressiveViaFormula(GA3E a, GA3E b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RegressiveViaFormula(r,s);
	cout << "u = RegressiveViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E LowerRightViaFormula(GA3E a, GA3E b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LowerRightViaFormula(r,s);
	cout << "u = LowerRightViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E Expander(const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Expander(r,s);
	cout << "u = Expander(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E Conserver(const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Conserver(r,s);
	cout << "u = Conserver(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E Shrinker(const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Shrinker(r,s);
	cout << "u = Shrinker(r,s) = " << u << "\n\n";

///////////////////////////////////////////////////////

// GA3E Symmetric(const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Symmetric(r,s);
	cout << "u = Symmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E AntiSymmetric(const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiSymmetric(r,s);
	cout << "u = AntiSymmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E Inner(const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Inner(r,s);
	cout << "u = Inner(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E LeftContraction (const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LeftContraction(r,s);
	cout << "u = LeftContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E RightContraction (const GA3E &a, const GA3E &b) ;

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	s = GA3E(b_q, b_x,b_y,b_z, b_xy,b_xz,b_yz, b_xyz); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RightContraction(r,s);
	cout << "u = RightContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// ex Determinant(GA3E A) ;

	X = GA3E(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	cout << "Determinant(X) = " << Determinant(X) << "\n\n";


//////////////////////////////////////////////////////

// GA3E Adjugate(GA3E V) ;

	X = GA3E(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Adjugate(X);
	cout << "u = Adjugate(X) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA3E Reciprocal(GA3E a) ;

	X = GA3E(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Reciprocal(X);
	cout << "u = Reciprocal(X) = " << u << "\n";
	s = X*u;
	cout << "X*(1/X) = " << s << "\n\n";

/////////////////////////////////  

	r = GA3E(a_q, a_x,a_y,a_z, a_xy,a_xz,a_yz, a_xyz);  
	det_ref = Determinant(r);
	int i;

	printf("Verify Magic transform preserves determinant\n");
	for(i=0; i<6; i++) {
		s = Magic(r,i);
		det = Determinant(s);
		if(det != det_ref) printf("Magic failure at index %d \n", i);
			else printf("+");
	}
	printf("\n\n");

/////////////////////////////////  

	printf("Verify Comp transform preserves determinant\n");
	for(i=0; i<32; i++) {
		s = Comp(r,i);
		det = Determinant(s);
		if(det != det_ref) printf("Comp failure at index %d \n", i);
			else printf("+");
	}
	printf("\n\n");

/////////////////////////////////  

	matrix N;

	GA3E MV, NV;

	MV = GA3E(3, 5,7,11, 13,17,19, 23);
	printf("MV = ");	PrintMV(MV);  printf("\n\n");

	N = GA3E_To_Matrix(MV);
	cout << "N = " << N << "\n";

	NV = Matrix_To_GA3E(N);
	printf("NV = ");	PrintMV(NV);  printf("\n\n");

/////////////////////////////////  

	return 0;
}




