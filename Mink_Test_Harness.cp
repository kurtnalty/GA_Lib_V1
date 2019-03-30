#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

#include "Mink_Routines.cp"

/////////////////////////////////// Demonstrate Mink Routines ////////////////////////////////


int main(void)
{
	
	Mink Zero;   Zero = Mink();

	symbol a_q("a");
	symbol a_x("b") ,  a_y("c") ,  a_z("d"),   a_t("e");
	symbol a_xy("f"),  a_xz("g"),  a_yz("h"),  a_xt("i"), a_yt("j"), a_zt("k");
	symbol a_xyz("l"), a_xyt("m"), a_xzt("n"), a_yzt("o");
	symbol a_xyzt("p");

	symbol b_q("A");
	symbol b_x("B") ,  b_y("C") ,  b_z("D"),   b_t("E");
	symbol b_xy("F"),  b_xz("G"),  b_yz("H"),  b_xt("I"), b_yt("J"), b_zt("K");
	symbol b_xyz("L"), b_xyt("M"), b_xzt("N"), b_yzt("O");
	symbol b_xyzt("P");

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

	Mink r,s,t,u,w;

	Mink X, Y, Z;

	X = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	Y = Mink( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	Z = Mink(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

	ex det_ref, det;
	det_ref = Determinant(r);

//////////////////////////////////////////////////////

// ostream &operator<<(ostream &ff, Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	cout << "\nr = " << r << "\n";

//////////////////////////////////////////////////////

// void PrintMV(Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	cout << "r = "; PrintMV(r); cout << "\n";

//////////////////////////////////////////////////////

// Mink OverBar(Mink a) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

	u = OverBar(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "OverBar(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*OverBar(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// Mink UnderBar(Mink a) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink Reverse(Mink w) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink Involution(Mink w) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink Transpose(Mink w) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink Conjugation(Mink w) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink CliffordConjugation(Mink w) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink Dual(Mink w) ; 

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink DorstDual(Mink a) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

	u = DorstDual(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "DorstDual(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*DorstDual(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// Mink DorstUnDual(Mink a) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  

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

// Mink operator+(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	u = r + s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r + s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// Mink operator-(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	u = r - s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r - s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator==(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r == s) cout << "r == s\n";   else cout << "r != s (as expected)\n";
	if (s == t) cout << "s == t (as expected)\n";   else cout << "s != t \n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator!=(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r != s) cout << "r != s (as expected)\n";   else cout << "r == s \n";
	if (s != t) cout << "s != t \n";   else cout << "s == t (as expected)\n";
	cout << "\n";

//////////////////////////////////////////////////////

Mink operator*(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r*s;
	cout << "u = r*s = " << u << "\n";

	u = r*s - s*r;
	if (u == Zero) cout << "Product is commutative\n"; else cout << "Product is non-commutative (as expected)\n";
	cout << "\n";

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
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

// Mink Product(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Product(r,s);
	cout << "u = Product(r,s) = " << u << "\n";

	u = Product(r,s) - Product(s,r);
	if (u == Zero) cout << "Product(r,s) is commutative\n"; else cout << "Product(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Product(Product(r,s), t) - Product(r,Product(s,t));
	if (u == Zero) cout << "Product is associative (as expected)\n"; else cout << "Product is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// Mink operator/(const Mink &u, const int i) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	cout << "r = "; PrintMV(r); cout << "\n";
	u = r/2;
	cout << "u = r/2 = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// Mink operator^(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r^s;
	cout << "u = r^s = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink Wedge(const Mink &u, const Mink &v) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Wedge(r,s);
	cout << "u = Wedge(r,s) = " << u << "\n";

	u = Wedge(r,s) - Wedge(s,r);
	if (u == Zero) cout << "Wedge(r,s) is commutative\n"; else cout << "Wedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Wedge(Wedge(r,s), t) - Wedge(r,Wedge(s,t));
	if (u == Zero) cout << "Wedge is associative (as expected)\n"; else cout << "Wedge is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// Mink AntiWedge(const Mink a, const Mink b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiWedge(r,s);
	cout << "u = AntiWedge(r,s) = " << u << "\n";

	u = AntiWedge(r,s) - AntiWedge(s,r);
	if (u == Zero) cout << "AntiWedge(r,s) is commutative\n"; else cout << "AntiWedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = AntiWedge(AntiWedge(r,s), t) - AntiWedge(r,AntiWedge(s,t));
	if (u == Zero) cout << "AntiWedge is associative (as expected)\n"; else cout << "AntiWedge is non-associative\n";
	cout << "\n";


//////////////////////////////////////////////////////

// Mink Regressive(Mink a, Mink b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Regressive(r,s);
	cout << "u = Regressive(r,s) = " << u << "\n";

	u = Regressive(r,s) - Regressive(s,r);
	if (u == Zero) cout << "Regressive(r,s) is commutative\n"; else cout << "Regressive(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Regressive(Regressive(r,s), t) - Regressive(r,Regressive(s,t));
	if (u == Zero) cout << "Regressive is associative (as expected)\n"; else cout << "Regressive is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// Mink RegressiveViaFormula(Mink a, Mink b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RegressiveViaFormula(r,s);
	cout << "u = RegressiveViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink LowerRightViaFormula(Mink a, Mink b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LowerRightViaFormula(r,s);
	cout << "u = LowerRightViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink Expander(const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LowerRightViaFormula(r,s);
	cout << "u = LowerRightViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink Conserver(const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Conserver(r,s);
	cout << "u = Conserver(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink Shrinker(const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Shrinker(r,s);
	cout << "u = Shrinker(r,s) = " << u << "\n\n";

///////////////////////////////////////////////////////

// Mink Symmetric(const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Symmetric(r,s);
	cout << "u = Symmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink AntiSymmetric(const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiSymmetric(r,s);
	cout << "u = AntiSymmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink Inner(const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Inner(r,s);
	cout << "u = Inner(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink LeftContraction (const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LeftContraction(r,s);
	cout << "u = LeftContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink RightContraction (const Mink &a, const Mink &b) ;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RightContraction(r,s);
	cout << "u = RightContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// ex Determinant(Mink A) ;

	X = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	cout << "Determinant(X) = " << Determinant(X) << "\n\n";


//////////////////////////////////////////////////////

// Mink Adjugate(Mink V) ;

	X = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Adjugate(X);
	cout << "u = Adjugate(X) = " << u << "\n\n";

//////////////////////////////////////////////////////

// Mink Reciprocal(Mink a) ;

	X = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Reciprocal(X);
	cout << "u = Reciprocal(X) = " << u << "\n";
	s = X*u;
	cout << "X*(1/X) = " << s << "\n\n";

/////////////////////////////////  

// Verify 64 Comps preserving the determinant

	int i;

	X = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	det_ref = Determinant(X);

	printf("Verifying 64 Determinant preserving Comp transforms\n");
	for (i=0; i<64; i++) {
		Y = Comp(X,i);
		det = Determinant(Y);
		if (det == det_ref) {
			printf("+");
		}
		else printf("Comp determinant mismatch at index %d \n", i);
	}
	printf("\n");


// Verify 72 Magic preserving the determinant

	X = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	det_ref = Determinant(X);

	printf("Verifying 72 Determinant preserving Magic transforms\n");
	for (i=0; i<72; i++) {
		Y = Magic(X,i);
		det = Determinant(Y);
		if (det == det_ref) {
			printf("+");
		}
		else printf("Magic determinant mismatch at index %d \n", i);
	}
	cout << "\n";

/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	X = Mink(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	Y = Mink( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	Z = Mink(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s = Mink(b_q,b_x,b_y,b_z,b_t,b_xy,b_xz,b_yz,b_xt,b_yt,b_zt,b_xyz,b_xyt,b_xzt,b_yzt,b_xyzt); 
	t = Mink(c_q,c_x,c_y,c_z,c_t,c_xy,c_xz,c_yz,c_xt,c_yt,c_zt,c_xyz,c_xyt,c_xzt,c_yzt,c_xyzt); 
	u = Mink(d_q,d_x,d_y,d_z,d_t,d_xy,d_xz,d_yz,d_xt,d_yt,d_zt,d_xyz,d_xyt,d_xzt,d_yzt,d_xyzt); 

*/





