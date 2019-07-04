// Routines for Geometric Algebra in Two Dimensional Euclidean spacetime
//
// Author: Kurt Nalty
// Release: 1.1
// Date: 27 October 2018
// License: Freeware
// Alternative License: BSD
//
//////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>

//////////////////////////////////////////////////////

typedef struct
{
	double q, x, y, xy;
} GA2E;

//////////////////////////////////////////////////////

GA2E GA2E_Zero(void)		// initializer
{
	GA2E a;

	a.q = 0.0; a.x = 0.0; a.y = 0.0; a.xy = 0.0;

	return a;

}

//////////////////////////////////////////////////////

GA2E GA2E_Set(double q, double x,  double y, double xy)	
{
	GA2E a;

	a.q = q;
	a.x = x; a.y = y;
	a.xy = xy;

	return a;

}

//////////////////////////////////////////////////////

// Necessary forward declarations

GA2E GA2E_Left_Contraction ( GA2E a,  GA2E b) ;
GA2E GA2E_Product( GA2E a,  GA2E b) ;


//////////////////////////////////////////////////////

void GA2E_PrintlnMV(GA2E v) 
{
	printf("(%10.3e,\n  %10.3e,%10.3e,\n  %10.3e ) \n",v.q, v.x,v.y, v.xy);
}

//////////////////////////////////////////////////////

void GA2E_PrintMV(GA2E v) 
{
	printf("(%10.3e,  %10.3e,%10.3e,  %10.3e ) ",v.q, v.x,v.y, v.xy);
}

//////////////////////////////////////////////////////

GA2E GA2E_OverBar(GA2E a)	 // in 3D, GA2E_OverBar and GA2E_UnderBar coincide
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

GA2E GA2E_UnderBar(GA2E a)   // in 3D, GA2E_OverBar and GA2E_UnderBar coincide
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

GA2E GA2E_Reverse(GA2E w)	
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

GA2E GA2E_Parity(GA2E w)	
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

GA2E GA2E_Transpose(GA2E w)	
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

GA2E GA2E_Conjugation(GA2E w)
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

GA2E GA2E_Clifford_Conjugation(GA2E w)
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

GA2E GA2E_Dual(GA2E w)   // return w*I_inv 
//     GA2E_Dual(r) = u = (d, c,-b, -a)

{
	GA2E v;
	GA2E I_inv;
	I_inv.xy = -1;
	v = GA2E_Product(w,I_inv);

	v.q =   w.xy;
	v.x =   w.y;
	v.y =  -w.x;
	v.xy = -w.q;

	return v;
}



//////////////////////////////////////////////////////

GA2E GA2E_Dorst_Dual(GA2E a)  // GA2E_Dorst_Dual(r) = u = (d, -c,b, -a)
{
	GA2E b;

//	GA2E I_inv;

//	I_inv = GA2E_Zero();
//	I_inv.xy = -1;

//	b = GA2E_Left_Contraction(a,I_inv);

	b.q =   a.xy;
	b.x =  -a.y;
	b.y =   a.x;
	b.xy = -a.q;

	return b;
}

//////////////////////////////////////////////////////

GA2E GA2E_Dorst_UnDual(GA2E a)  //GA2E_Dorst_UnDual(r) = u = (-d, c,-b, a)
{
	GA2E b;
//	GA2E I;	I.xy = 1;

//	b = GA2E_Left_Contraction(a,I);

	b.q =  -a.xy;
	b.x =   a.y;
	b.y =  -a.x;
	b.xy =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA2E GA2E_Add( GA2E u,  GA2E v)
{
	GA2E w;
	w.q =   u.q   + v.q  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.xy  = u.xy  + v.xy ;

	return w;
}

//////////////////////////////////////////////////////

GA2E GA2E_Subtract( GA2E u,  GA2E v)
{
	GA2E w;
	w.q =   u.q   - v.q  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.xy  = u.xy  - v.xy ;

	return w;
}

//////////////////////////////////////////////////////

int GA2E_Equal( GA2E u,  GA2E v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.x ==v.x )&&
			(u.y ==v.y )&&

			(u.xy==v.xy) ;
	return result;
}

//////////////////////////////////////////////////////

int GA2E_Not_Equal( GA2E u,  GA2E v)
{
	int result;
	result = 	(u.q !=v.q )||

			(u.x !=v.x )||
			(u.y !=v.y )||

			(u.xy!=v.xy) ;
	return result;
}

//////////////////////////////////////////////////////

GA2E GA2E_Product( GA2E a,  GA2E b)
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

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Multiply_By_Constant( GA2E u, double a)
{
	GA2E w;

	w.q =   u.q*a;

	w.x   = u.x*a;
	w.y   = u.y*a;

	w.xy  = u.xy*a;

	return w;
}

//////////////////////////////////////////////////////

GA2E GA2E_Divide_By_Constant( GA2E u, double a)
{
	GA2E w;

	w.q =   u.q/a;

	w.x   = u.x/a;
	w.y   = u.y/a;

	w.xy  = u.xy/a;

	return w;
}

//////////////////////////////////////////////////////

GA2E GA2E_Wedge( GA2E a,  GA2E b) {

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

	return c;

}

//////////////////////////////////////////////////////

GA2E GA2E_AntiWedge( GA2E a,  GA2E b)
{
	GA2E c;
/*

Lengyel's GA2E_AntiWedge
 V    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     q    
 x    |  0     0     q     x    
 y    |  0    -q     0     y    
 xy   |  q     x     y     xy   

*/

//	c = GA2E_OverBar(GA2E_Wedge(GA2E_UnderBar(a),GA2E_UnderBar(b)));

	c.q   =  + a.q *b.xy  + a.x *b.y - a.y*b.x + a.xy*b.q   ; 
	c.x   =  + a.x *b.xy  + a.xy*b.x   ; 
	c.y   =  + a.y *b.xy  + a.xy*b.y   ; 
	c.xy  =  + a.xy*b.xy  ; 

	return c;
}

//////////////////////////////////////////////////////


GA2E GA2E_Regressive(GA2E a, GA2E b) 
{
	GA2E c;	

/*

Hestenes' GA2E_Regressive (differs from GA2E_AntiWedge in signs)
 V    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     q    
 x    |  0     0    -q     x    
 y    |  0     q     0     y    
 xy   |  q     x     y     xy   

*/

	c.q = a.q*b.xy - a.x*b.y + a.xy*b.q + a.y*b.x ;
	c.x =  + (a.x*b.xy + a.xy*b.x) ;
	c.y =  + (a.xy*b.y + a.y*b.xy) ;
	c.xy = + a.xy*b.xy ;

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Regressive_Via_Formula(GA2E a, GA2E b)
{
	GA2E c;

	GA2E Ixy,I_inv;
	
	Ixy = GA2E_Zero();	I_inv = GA2E_Zero();
	Ixy.xy = 1;	I_inv.xy = -1;

	c = GA2E_Product(GA2E_Wedge(GA2E_Product(a,I_inv),GA2E_Product(b,I_inv)),Ixy);

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Lower_Right_Via_Formula(GA2E a, GA2E b)  // Like a mirror of GA2E_Wedge along rising diagonal
{

/*

GA2E_Lower_Right_Via_Formula(Blade[i],Blade[j])
 ?    |  q     x     y     xy   
------------------------------
 q    |  0     0     0     xy   
 x    |  0     0     xy    y    
 y    |  0    -xy    0    -x    
 xy   |  xy   -y     x    -q    

*/
// GA2E_Wedge(a*I_inv,I*b);

	GA2E c;

	GA2E Ixy,I_inv;
	
	Ixy = GA2E_Zero();	I_inv = GA2E_Zero();
	Ixy.xy = 1;	I_inv.xy = -1;

	c = GA2E_Wedge(GA2E_Product(a,I_inv),GA2E_Product(Ixy,b));

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Expander( GA2E a,  GA2E b)
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

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Conserver( GA2E a,  GA2E b)
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

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Shrinker( GA2E a,  GA2E b)
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

// GA2E_Shrinker equation set KN
	c.q    =  + a.x*b.x + a.y*b.y - a.xy*b.xy; 
	c.x    =  - a.y*b.xy + a.xy*b.y ; 
	c.y    =  + a.x*b.xy - a.xy*b.x ; 
	c.xy   =  0 ; 

	return c;
}

///////////////////////////////////////////////////////

GA2E GA2E_Symmetric( GA2E a,  GA2E b)
{ 
/*

GA2E_Symmetric GA2E_Product
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

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_AntiSymmetric( GA2E a,  GA2E b)
{ 
/*

GA2E_AntiSymmetric GA2E_Product
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

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Inner( GA2E a,  GA2E b)
{ 
/*

GA2E_Inner GA2E_Product
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

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Left_Contraction ( GA2E a,  GA2E b)
{ 
/*

GA2E_Left_Contraction GA2E_Product
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

	return c;
}

//////////////////////////////////////////////////////

GA2E GA2E_Right_Contraction ( GA2E a,  GA2E b)
{ 
/*

GA2E_Right_Contraction GA2E_Product
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

	return c;
}

//////////////////////////////////////////////////////

double GA2E_Determinant(GA2E A) {

	double e, det;
/*
	r = (a, b,c, d)
	s = GA2E_UnderBar(r) = (d, c,-b, a)
	t = r*GA2E_UnderBar(r) = (0, 0,0, -c^2+a^2+d^2-b^2)

*/
	e = A.q*A.q - A.x*A.x - A.y*A.y + A.xy*A.xy;
	det = e*e;

	return(det);
}

//////////////////////////////////////////////////////

GA2E GA2E_Adjugate(GA2E a)
{

	GA2E u;

//GA2E_Product(GA2E_Reverse(r),GA2E_Conjugation(GA2E_Product(r,GA2E_Reverse(r)))) = GA2E_Adjugate(r); KN

	u = GA2E_Product(GA2E_Reverse(a),GA2E_Conjugation(GA2E_Product(a,GA2E_Reverse(a))));

	return u;

}

//////////////////////////////////////////////////////

GA2E GA2E_Reciprocal(GA2E a) 
{
	double b;
	GA2E c,d;

	b = GA2E_Determinant(a);

	c = GA2E_Adjugate(a);

	d.q = (c.q/b);

	d.x = (c.x/b);
	d.y = (c.y/b);

	d.xy = (c.xy/b);

	return d;
}

//////////////////////////////////////////////////////

// int	main(void) { return 0; }





