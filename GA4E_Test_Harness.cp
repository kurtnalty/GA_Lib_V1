#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

#include "GA4E_Routines.cp"

/////////////////////////////////// Demonstrate GA4E Routines ////////////////////////////////


int main(void)
{
	
	GA4E Zero;   Zero = GA4E();

	symbol a_q("a");
	symbol a_x("b") ,  a_y("c") ,  a_z("d"),   a_w("e");
	symbol a_xy("f"),  a_xz("g"),  a_yz("h"),  a_xw("i"), a_yw("j"), a_zw("k");
	symbol a_xyz("l"), a_xyw("m"), a_xzw("n"), a_yzw("o");
	symbol a_xyzw("p");

	symbol b_q("A");
	symbol b_x("B") ,  b_y("C") ,  b_z("D"),   b_w("E");
	symbol b_xy("F"),  b_xz("G"),  b_yz("H"),  b_xw("I"), b_yw("J"), b_zw("K");
	symbol b_xyz("L"), b_xyw("M"), b_xzw("N"), b_yzw("O");
	symbol b_xyzw("P");

	symbol c_q("c.a");
	symbol c_x("c.b") ,  c_y("c.c") ,  c_z("c.d"),   c_w("c.e");
	symbol c_xy("c.f"),  c_xz("c.g"),  c_yz("c.h"),  c_xw("c.i"), c_yw("c.j"), c_zw("c.k");
	symbol c_xyz("c.l"), c_xyw("c.m"), c_xzw("c.n"), c_yzw("c.o");
	symbol c_xyzw("c.p");

	symbol d_q("d.a");
	symbol d_x("d.b") ,  d_y("d.c") ,  d_z("d.d"),   d_w("d.e");
	symbol d_xy("d.f"),  d_xz("d.g"),  d_yz("d.h"),  d_xw("d.i"), d_yw("d.j"), d_zw("d.k");
	symbol d_xyz("d.l"), d_xyw("d.m"), d_xzw("d.n"), d_yzw("d.o");
	symbol d_xyzw("d.p");

	GA4E r,s,t,u,w;

	GA4E X, Y, Z;

	X = GA4E(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	Y = GA4E( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	Z = GA4E(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 
	u = GA4E(d_q,d_x,d_y,d_z,d_w,d_xy,d_xz,d_yz,d_xw,d_yw,d_zw,d_xyz,d_xyw,d_xzw,d_yzw,d_xyzw); 

	ex det_ref, det;
	det_ref = Determinant(r);

//////////////////////////////////////////////////////

// ostream &operator<<(ostream &ff, GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	cout << "\nr = " << r << "\n";

//////////////////////////////////////////////////////

// void PrintMV(GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	cout << "r = "; PrintMV(r); cout << "\n";

//////////////////////////////////////////////////////

// GA4E OverBar(GA4E a) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

	u = OverBar(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "OverBar(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*OverBar(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA4E UnderBar(GA4E a) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E Reverse(GA4E w) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E Involution(GA4E w) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E Transpose(GA4E w) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E Conjugation(GA4E w) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E CliffordConjugation(GA4E w) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E Dual(GA4E w) ; 

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E DorstDual(GA4E a) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

	u = DorstDual(r);

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "DorstDual(r) = "; PrintMV(u); cout << "\n";
	det = Determinant(u);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	s = r*u;
	cout << "r*DorstDual(r) = " << s  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA4E DorstUnDual(GA4E a) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  

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

// GA4E operator+(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	u = r + s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r + s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA4E operator-(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	u = r - s;

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "r - s = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator==(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r == s) cout << "r == s\n";   else cout << "r != s (as expected)\n";
	if (s == t) cout << "s == t (as expected)\n";   else cout << "s != t \n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator!=(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 

	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	if (r != s) cout << "r != s (as expected)\n";   else cout << "r == s \n";
	if (s != t) cout << "s != t \n";   else cout << "s == t (as expected)\n";
	cout << "\n";

//////////////////////////////////////////////////////

GA4E operator*(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r*s;
	cout << "u = r*s = " << u << "\n";

	u = r*s - s*r;
	if (u == Zero) cout << "Product is commutative\n"; else cout << "Product is non-commutative (as expected)\n";
	cout << "\n";

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 
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

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 


//////////////////////////////////////////////////////

// GA4E Product(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Product(r,s);
	cout << "u = Product(r,s) = " << u << "\n";

	u = Product(r,s) - Product(s,r);
	if (u == Zero) cout << "Product(r,s) is commutative\n"; else cout << "Product(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Product(Product(r,s), t) - Product(r,Product(s,t));
	if (u == Zero) cout << "Product is associative (as expected)\n"; else cout << "Product is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA4E operator/(const GA4E &u, const int i) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	cout << "r = "; PrintMV(r); cout << "\n";
	u = r/2;
	cout << "u = r/2 = "; PrintMV(u); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA4E operator^(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = r^s;
	cout << "u = r^s = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E Wedge(const GA4E &u, const GA4E &v) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Wedge(r,s);
	cout << "u = Wedge(r,s) = " << u << "\n";

	u = Wedge(r,s) - Wedge(s,r);
	if (u == Zero) cout << "Wedge(r,s) is commutative\n"; else cout << "Wedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Wedge(Wedge(r,s), t) - Wedge(r,Wedge(s,t));
	if (u == Zero) cout << "Wedge is associative (as expected)\n"; else cout << "Wedge is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA4E AntiWedge(const GA4E a, const GA4E b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiWedge(r,s);
	cout << "u = AntiWedge(r,s) = " << u << "\n";

	u = AntiWedge(r,s) - AntiWedge(s,r);
	if (u == Zero) cout << "AntiWedge(r,s) is commutative\n"; else cout << "AntiWedge(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = AntiWedge(AntiWedge(r,s), t) - AntiWedge(r,AntiWedge(s,t));
	if (u == Zero) cout << "AntiWedge is associative (as expected)\n"; else cout << "AntiWedge is non-associative\n";
	cout << "\n";


//////////////////////////////////////////////////////

// GA4E Regressive(GA4E a, GA4E b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Regressive(r,s);
	cout << "u = Regressive(r,s) = " << u << "\n";

	u = Regressive(r,s) - Regressive(s,r);
	if (u == Zero) cout << "Regressive(r,s) is commutative\n"; else cout << "Regressive(r,s) is non-commutative (as expected)\n";
	cout << "\n";

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	cout << "t = "; PrintMV(t); cout << "\n";
	u = Regressive(Regressive(r,s), t) - Regressive(r,Regressive(s,t));
	if (u == Zero) cout << "Regressive is associative (as expected)\n"; else cout << "Regressive is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA4E RegressiveViaFormula(GA4E a, GA4E b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RegressiveViaFormula(r,s);
	cout << "u = RegressiveViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E LowerRightViaFormula(GA4E a, GA4E b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LowerRightViaFormula(r,s);
	cout << "u = LowerRightViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E Expander(const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LowerRightViaFormula(r,s);
	cout << "u = LowerRightViaFormula(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E Conserver(const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Conserver(r,s);
	cout << "u = Conserver(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E Shrinker(const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Shrinker(r,s);
	cout << "u = Shrinker(r,s) = " << u << "\n\n";

///////////////////////////////////////////////////////

// GA4E Symmetric(const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Symmetric(r,s);
	cout << "u = Symmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E AntiSymmetric(const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = AntiSymmetric(r,s);
	cout << "u = AntiSymmetric(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E Inner(const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = Inner(r,s);
	cout << "u = Inner(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E LeftContraction (const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = LeftContraction(r,s);
	cout << "u = LeftContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E RightContraction (const GA4E &a, const GA4E &b) ;

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	cout << "r = "; PrintMV(r); cout << "\n";
	cout << "s = "; PrintMV(s); cout << "\n";
	u = RightContraction(r,s);
	cout << "u = RightContraction(r,s) = " << u << "\n\n";

//////////////////////////////////////////////////////

// ex Determinant(GA4E A) ;

	X = GA4E(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	cout << "Determinant(X) = " << Determinant(X) << "\n\n";


//////////////////////////////////////////////////////

// GA4E Adjugate(GA4E V) ;

	X = GA4E(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Adjugate(X);
	cout << "u = Adjugate(X) = " << u << "\n\n";

//////////////////////////////////////////////////////

// GA4E Reciprocal(GA4E a) ;

	X = GA4E(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	cout << "X = "; PrintMV(X); cout << "\n";
	u = Reciprocal(X);
	cout << "u = Reciprocal(X) = " << u << "\n";
	s = X*u;
	cout << "X*(1/X) = " << s << "\n\n";

/////////////////////////////////  

// Verify 64 Comps preserving the determinant

	int i;

	X = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
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


// Verify 120 Magic preserving the determinant

	X = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	det_ref = Determinant(X);

	printf("Verifying 120 Determinant preserving Magic transforms\n");
	for (i=0; i<120; i++) {
		Y = Magic(X,i);
		det = Determinant(Y);
		if (det == det_ref) {
			printf("+");
		}
		else printf("Magic determinant mismatch at index %d \n", i);
	}
	printf("\n");



/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	X = GA4E(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	Y = GA4E( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	Z = GA4E(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

	r = GA4E(a_q,a_x,a_y,a_z,a_w,a_xy,a_xz,a_yz,a_xw,a_yw,a_zw,a_xyz,a_xyw,a_xzw,a_yzw,a_xyzw);  
	s = GA4E(b_q,b_x,b_y,b_z,b_w,b_xy,b_xz,b_yz,b_xw,b_yw,b_zw,b_xyz,b_xyw,b_xzw,b_yzw,b_xyzw); 
	t = GA4E(c_q,c_x,c_y,c_z,c_w,c_xy,c_xz,c_yz,c_xw,c_yw,c_zw,c_xyz,c_xyw,c_xzw,c_yzw,c_xyzw); 
	u = GA4E(d_q,d_x,d_y,d_z,d_w,d_xy,d_xz,d_yz,d_xw,d_yw,d_zw,d_xyz,d_xyw,d_xzw,d_yzw,d_xyzw); 

*/





