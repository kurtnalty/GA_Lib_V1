// Routines for Geometric Algebra in Three Dimensional Euclidean spacetime
//
// Author: Kurt Nalty
// Release: 1.0
// Date: 12 March 2018
// License: Freeware
// Alternative License: BSD
//
//////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

//////////////////////////////////////////////////////

struct GA2E{	// 0 . 1 2 . 3
	ex q,  x,y, xy;
	GA2E() {q = 0; x = 0; y = 0; xy = 0;}
	GA2E(ex qq, ex xx, ex yy, ex xxyy) 
		{q = qq; x = xx; y = yy; xy = xxyy;}
};


//////////////////////////////////////////////////////

// Necessary forward declarations

GA2E LeftContraction (const GA2E &a, const GA2E &b) ;
GA2E Product(const GA2E &a, const GA2E &b) ;

//////////////////////////////////////////////////////

ostream &operator<<(ostream &ff, GA2E &v) 
{
	return ff << "\n(" 
		<< v.q << ", "
		<< v.x << "," 
		<< v.y << ", "  
		<< v.xy << ")\n";
}


//////////////////////////////////////////////////////

void PrintMV(GA2E &v) 
{
	cout 	<< "( "
		<< v.q << ", " 
		<< v.x << "," 
		<< v.y << ", " 
		<< v.xy << ")";
}

//////////////////////////////////////////////////////

GA2E OverBar(GA2E a)	 // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA2E b;		// a blade wedge b blade = pseudovector blade xy

	b.q     =  a.xy;	// xy
	b.x     = -a.y;
	b.y     =  a.x;
	b.xy   =   a.q;

	return b;
}

//////////////////////////////////////////////////////

GA2E UnderBar(GA2E a)   // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA2E b;		// b blade wedge a blade = pseudovector blade

	b.q     =  a.xy;	// xy
	b.x     =  a.y;
	b.y     = -a.x;
	b.xy   =   a.q;

	return b;
}

//////////////////////////////////////////////////////

GA2E Reverse(GA2E w)	
// (A,  B, C,-D) 
{
	GA2E v;
	v.q =  w.q;
	v.x = w.x;
	v.y = w.y;
	v.xy = -w.xy;

	return v;
}

//////////////////////////////////////////////////////

GA2E Involution(GA2E w)	
// Corresponds to P parity transform: x -> -x, y -> -y,
// (A, -B,-C, D) 
{
	GA2E v;
	v.q =  w.q;
	v.x = -w.x;
	v.y = -w.y;
	v.xy = w.xy;

	return v;
}

//////////////////////////////////////////////////////

GA2E Transpose(GA2E w)	
//(A,  B, C,-D)  
{
	GA2E v;

	v.q =  w.q;
	v.x =  w.x;
	v.y =  w.y;
	v.xy = -w.xy;

	return v;
}

//////////////////////////////////////////////////////

GA2E Conjugation(GA2E w)
// (A, -B,-C,-D)
{
	GA2E v;
	v.q =  w.q;
	v.x = -w.x;
	v.y = -w.y;
	v.xy = -w.xy;

	return v;
}

//////////////////////////////////////////////////////

GA2E CliffordConjugation(GA2E w)
// (A, -B,-C,-D, -E,-F,-G,  H) 
{
	GA2E v;
	v.q =  w.q;
	v.x = -w.x;
	v.y = -w.y;
	v.xy = -w.xy;

	return v;
}

//////////////////////////////////////////////////////

GA2E Dual(GA2E w)   // return w*I_inv 
//     Dual(r) = u = (d, c,-b, -a)

{
	GA2E v;
	GA2E I_inv;
	I_inv.xy = -1;
	v = Product(w,I_inv);

	v.q =   w.xy;
	v.x =   w.y;
	v.y =  -w.x;
	v.xy = -w.q;

	return v;
}



//////////////////////////////////////////////////////

GA2E DorstDual(GA2E a)  // DorstDual(r) = u = (d, -c,b, -a)
{
	GA2E b;
	GA2E I_inv;	I_inv.xy = -1;

	b = LeftContraction(a,I_inv);

	b.q =   a.xy;
	b.x =  -a.y;
	b.y =   a.x;
	b.xy = -a.q;

	return b;
}

//////////////////////////////////////////////////////

GA2E DorstUnDual(GA2E a)  //DorstUnDual(r) = u = (-d, c,-b, a)
{
	GA2E b;
	GA2E I;	I.xy = 1;

	b = LeftContraction(a,I);

	b.q =  -a.xy;
	b.x =   a.y;
	b.y =  -a.x;
	b.xy =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA2E operator+(const GA2E &u, const GA2E &v)
{
	GA2E w;
	w.q =   u.q   + v.q  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.xy  = u.xy  + v.xy ;

	return w;
}

//////////////////////////////////////////////////////

GA2E operator-(const GA2E &u, const GA2E &v)
{
	GA2E w;
	w.q =   u.q   - v.q  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.xy  = u.xy  - v.xy ;

	return w;
}

//////////////////////////////////////////////////////

int operator==(const GA2E &u, const GA2E &v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.x ==v.x )&&
			(u.y ==v.y )&&

			(u.xy==v.xy) ;
	return result;
}

//////////////////////////////////////////////////////

int operator!=(const GA2E &u, const GA2E &v)
{
	int result;
	result = 	(u.q !=v.q )||

			(u.x !=v.x )||
			(u.y !=v.y )||

			(u.xy!=v.xy) ;
	return result;
}

//////////////////////////////////////////////////////

GA2E operator*(const GA2E &a, const GA2E &b) {
	GA2E c;

/*

 *    |  q     x     y     xy
-----------------------------
 q    |  q     x     y     xy 
 x    |  x     q     xy    y    
 y    |  y    -xy    q    -x    
 xy   |  xy   -y     x    -q  

*/
	c.q   =  + a.q*b.q   + a.x*b.x   + a.y*b.y  - a.xy*b.xy ;
	c.x   =  + a.q*b.x   + a.x*b.q   - a.y*b.xy + a.xy*b.y  ;
	c.y   =  + a.q*b.y   + a.x*b.xy  + a.y*b.q  - a.xy*b.x  ;
	c.xy  =  + a.q*b.xy  + a.x*b.y   - a.y*b.x  + a.xy*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E Product(const GA2E &a, const GA2E &b)
{

	GA2E c;

/*

 *    |  q     x     y     xy
-----------------------------
 q    |  q     x     y     xy 
 x    |  x     q     xy    y    
 y    |  y    -xy    q    -x    
 xy   |  xy   -y     x    -q  

*/
	c.q   =  + a.q*b.q   + a.x*b.x   + a.y*b.y  - a.xy*b.xy ;
	c.x   =  + a.q*b.x   + a.x*b.q   - a.y*b.xy + a.xy*b.y  ;
	c.y   =  + a.q*b.y   + a.x*b.xy  + a.y*b.q  - a.xy*b.x  ;
	c.xy  =  + a.q*b.xy  + a.x*b.y   - a.y*b.x  + a.xy*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E operator*(const GA2E &v, const ex &u) {
	GA2E w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.xy   = expand(u*v.xy  ) ;

	return w;
}

//////////////////////////////////////////////////////

GA2E operator*(const ex &u, const GA2E &v) {
	GA2E w;

	w.q    = expand(u*v.q   ) ;
	w.x    = expand(u*v.x   ) ;
	w.y    = expand(u*v.y   ) ;
	w.xy   = expand(u*v.xy  ) ;

	return w;
}


//////////////////////////////////////////////////////

GA2E operator/(const GA2E &u, const int i)
{
	GA2E w;

	w.q =   u.q/i;

	w.x   = u.x/i;
	w.y   = u.y/i;

	w.xy  = u.xy/i;

	return w;
}

//////////////////////////////////////////////////////

GA2E operator^(const GA2E &a, const GA2E &b) {

	GA2E c;

/*

 ^    |  q     x     y     xy   
------------------------------
 q    |  q     x     y     xy   
 x    |  x     0     xy    0    
 y    |  y    -xy    0     0    
 xy   |  xy    0     0     0    

*/

	c.q   =  + a.q*b.q  ;
	c.x   =  + a.q*b.x   + a.x*b.q  ;
	c.y   =  + a.q*b.y   + a.y*b.q  ;
	c.xy  =  + a.q*b.xy  + a.x*b.y  - a.y*b.x  + a.xy*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E Wedge(const GA2E &a, const GA2E &b) {

	GA2E c;


/*

 ^    |  q     x     y     xy
-----------------------------
 q    |  q     x     y     xy  
 x    |  x     0     xy    0     
 y    |  y    -xy    0     0   
 xy   |  xy    0     0     0     

*/

	c.q   =  + a.q*b.q  ;
	c.x   =  + a.q*b.x   + a.x*b.q  ;
	c.y   =  + a.q*b.y   + a.y*b.q  ;
	c.xy  =  + a.q*b.xy  + a.x*b.y  - a.y*b.x  + a.xy*b.q  ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;

}

//////////////////////////////////////////////////////

GA2E AntiWedge(const GA2E a, const GA2E b)
{
	GA2E c;
/*

Lengyel's AntiWedge
 V    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     q    
 x    |  0     0     q     x    
 y    |  0    -q     0     y    
 xy   |  q     x     y     xy   

*/

//	c = OverBar(Wedge(UnderBar(a),UnderBar(b)));

	c.q   =  + a.q *b.xy  + a.x *b.y - a.y*b.x + a.xy*b.q   ; 
	c.x   =  + a.x *b.xy  + a.xy*b.x   ; 
	c.y   =  + a.y *b.xy  + a.xy*b.y   ; 
	c.xy  =  + a.xy*b.xy  ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////


GA2E Regressive(GA2E a, GA2E b) 
{
	GA2E c;	

/*

Hestenes' Regressive (differs from AntiWedge in signs)
 V    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     q    
 x    |  0     0    -q     x    
 y    |  0     q     0     y    
 xy   |  q     x     y     xy   

*/

	GA2E I, I_inv;
	I.xy = 1;
	I_inv.xy = -1;

//	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = a.q*b.xy - a.x*b.y + a.xy*b.q + a.y*b.x ;
	c.x =  + (a.x*b.xy + a.xy*b.x) ;
	c.y =  + (a.xy*b.y + a.y*b.xy) ;
	c.xy = + a.xy*b.xy ;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);

	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E RegressiveViaFormula(GA2E a, GA2E b)
{
	GA2E c;

	GA2E I,I_inv;
	
	I.xy = 1;	I_inv.xy = -1;

	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = expand(c.q);

	c.x = expand(c.x);
	c.y = expand(c.y);

	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E LowerRightViaFormula(GA2E a, GA2E b)  // Like a mirror of Wedge along rising diagonal
{

/*

LowerRightViaFormula(Blade[i],Blade[j])
 ?    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     xy   
 x    |  0     0     xy    y    
 y    |  0    -xy    0    -x    
 xy   |  xy   -y     x    -q    

*/
// Wedge(a*I_inv,I*b);

	GA2E c;

	GA2E I,I_inv,d,e,f;
	
	I.xy = 1;	I_inv.xy = -1;

	c = Wedge(a*I_inv,I*b);
	

	return c;
}

//////////////////////////////////////////////////////

GA2E Expander(const GA2E &a, const GA2E &b)
{ 
/*

Terms with increased rank
 >    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     0    
 x    |  0     0     xy    0    
 y    |  0    -xy    0     0    
 xy   |  0     0     0     0    

*/
	GA2E c;

	c.q    =  0 ; 
	c.x    =  0 ; 
	c.y    =  0 ; 
	c.xy   =  + a.x*b.y - a.y*b.x    ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E Conserver(const GA2E &a, const GA2E &b)
{  

/*

Terms with preserved rank
 =    |  q     x     y     xy   
------------------------------
 q    |  q     x     y     xy   
 x    |  x     0     0     0    
 y    |  y     0     0     0    
 xy   |  xy    0     0     0    

*/
	GA2E c;

	c.q    =  + a.q*b.q    ; 
	c.x    =  + a.q*b.x    + a.x*b.q    ; 
	c.y    =  + a.q*b.y    + a.y*b.q    ; 	
	c.xy   =  + a.q*b.xy   + a.xy*b.q ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E Shrinker(const GA2E &a, const GA2E &b)
{

/*

Terms with reduced rank
 <    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     0    
 x    |  0     q     0     y    
 y    |  0     0     q    -x    
 xy   |  0    -y     x    -q    

*/
	GA2E c;

// Shrinker equation set KN
	c.q    =  + a.x*b.x + a.y*b.y - a.xy*b.xy; 
	c.x    =  - a.y*b.xy + a.xy*b.y ; 
	c.y    =  + a.x*b.xy - a.xy*b.x ; 
	c.xy   =  0 ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

///////////////////////////////////////////////////////

GA2E Symmetric(const GA2E &a, const GA2E &b)
{ 
/*

Symmetric Product
 ?    |  q     x     y     xy   
------------------------------
 q    |  q     x     y     xy   
 x    |  x     q     0     0    
 y    |  y     0     q     0    
 xy   |  xy    0     0    -q    

*/
	GA2E c;

//	c = (a*b + b*a)/2;

	c.q   =  + a.q*b.q  + a.x *b.x + a.y*b.y - a.xy*b.xy  ; 
	c.x   =  + a.q*b.x  + a.x *b.q   ; 
	c.y   =  + a.q*b.y  + a.y *b.q   ; 
	c.xy  =  + a.q*b.xy + a.xy*b.q   ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E AntiSymmetric(const GA2E &a, const GA2E &b)
{ 
/*

AntiSymmetric Product
 ?    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     0    
 x    |  0     0     xy    y    
 y    |  0    -xy    0    -x    
 xy   |  0    -y     x     0    

*/
	GA2E c;

//	c = (a*b - b*a)/2;

	c.q   = 0;
	c.x   =  - a.y*b.xy + a.xy*b.y   ; 
	c.y   =  + a.x*b.xy - a.xy*b.x   ; 
	c.xy  =  + a.x*b.y  - a.y *b.x   ; 

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E Inner(const GA2E &a, const GA2E &b)
{ 
/*

Inner Product
 ?    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     0    
 x    |  0     q     0     y    
 y    |  0     0     q    -x    
 xy   |  0    -y     x    -q    

*/
	GA2E c;

	c.q = a.x*b.x - a.xy*b.xy + a.y*b.y ;
	c.x = + (a.xy*b.y - a.y*b.xy) ;
	c.y = + (a.x*b.xy - a.xy*b.x) ;
	c.xy = 0 ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E LeftContraction (const GA2E &a, const GA2E &b)
{ 
/*

LeftContraction Product
 ?    |  q     x     y     xy   
------------------------------
 q    |  q     x     y     xy   
 x    |  0     q     0     y    
 y    |  0     0     q    -x    
 xy   |  0     0     0    -q    

*/
	GA2E c;

	c.q = a.q*b.q + a.x*b.x - a.xy*b.xy + a.y*b.y ;
	c.x = + (a.q*b.x - a.y*b.xy) ;
	c.y = + (a.q*b.y + a.x*b.xy) ;
	c.xy = + (a.q*b.xy) ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

GA2E RightContraction (const GA2E &a, const GA2E &b)
{ 
/*

RightContraction Product
 ?    |  q     x     y     xy   
------------------------------
 q    |  q     0     0     0    
 x    |  x     q     0     0    
 y    |  y     0     q     0    
 xy   |  xy   -y     x    -q    

*/
	GA2E c;

	c.q = a.q*b.q + a.x*b.x - a.xy*b.xy + a.y*b.y ;
	c.x = + ( a.x*b.q + a.xy*b.y) ;
	c.y = + (-a.xy*b.x + a.y*b.q) ;
	c.xy = + ( a.xy*b.q) ;

	c.q = expand(c.q);
	c.x = expand(c.x);
	c.y = expand(c.y);
	c.xy = expand(c.xy);

	return c;
}

//////////////////////////////////////////////////////

ex Determinant(GA2E A) {

	ex e, det;
/*
	r = (a, b,c, d)
	s = UnderBar(r) = (d, c,-b, a)
	t = r*UnderBar(r) = (0, 0,0, -c^2+a^2+d^2-b^2)
*/
	e = A.q*A.q - A.x*A.x - A.y*A.y + A.xy*A.xy;
	det = expand(e*e);

	return(det);
}

//////////////////////////////////////////////////////

GA2E Adjugate(GA2E a)
{

	GA2E u;

//Product(Reverse(r),Conjugation(Product(r,Reverse(r)))) = Adjugate(r); KN

	u = Product(Reverse(a),Conjugation(Product(a,Reverse(a))));

	u.q = expand(u.q);

	u.x = expand(u.x);
	u.y = expand(u.y);

	u.xy = expand(u.xy);


	return u;

}

//////////////////////////////////////////////////////

GA2E Reciprocal(GA2E a) // KN
{
	ex b;
	GA2E c,d;

	b = Determinant(a);

	c = Adjugate(a);

	d.q = expand(c.q/b);

	d.x = expand(c.x/b);
	d.y = expand(c.y/b);

	d.xy = expand(c.xy/b);

	return d;
}

//////////////////////////////////////////////////////

// int	main(void) { return 0; }





