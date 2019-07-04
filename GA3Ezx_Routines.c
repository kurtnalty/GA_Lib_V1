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
	double q;
	double x,  y,  z;
	double xy, zx, yz;
	double xyz;
} GA3Ezx;


//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Zero(void)		// initializer
{
	GA3Ezx a;

	a.q = 0.0;
	a.x = 0.0; a.y = 0.0; a.z = 0.0;
	a.xy = 0.0; a.zx = 0.0; a.yz = 0.0;
	a.xyz = 0.0;

	return a;

}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Set(	double q,
	double x,  double y,  double z,
	double xy, double zx, double yz,
	double xyz)		// initializer
{
	GA3Ezx a;

	a.q = q;
	a.x = x; a.y = y; a.z = z;
	a.xy = xy; a.zx = zx; a.yz = yz;
	a.xyz = xyz;

	return a;

}

//////////////////////////////////////////////////////

// Necessary forward declarations

GA3Ezx GA3Ezx_LeftContraction (GA3Ezx a, GA3Ezx b) ;
GA3Ezx GA3Ezx_Product(GA3Ezx a, GA3Ezx b) ;

//////////////////////////////////////////////////////

void GA3Ezx_PrintlnMV(GA3Ezx v) 
{
	printf("(%10.3e,\n  %10.3e,%10.3e,%10.3e,\n  %10.3e,%10.3e,%10.3e,\n  %10.3e) \n",v.q, v.x,v.y,v.z, v.xy,v.zx,v.yz, v.xyz);
}

//////////////////////////////////////////////////////

void GA3Ezx_PrintMV(GA3Ezx v) 
{
//	printf("(%e,  %e,%e,%e,  %e,%e,%e,  %e) ",v.q, v.x,v.y,v.z, v.xy,v.zx,v.yz, v.xyz);
	printf("(%10.3e,  %10.3e,%10.3e,%10.3e,  %10.3e,%10.3e,%10.3e,  %10.3e) ",v.q, v.x,v.y,v.z, v.xy,v.zx,v.yz, v.xyz);
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_OverBar(GA3Ezx a)	 // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA3Ezx b;		// a blade wedge b blade = pseudovector blade
// OverBar(Blade[i]) = Reciprocal(Blade[i])*xyz ;

	b.q     =  a.xyz;	// xyz

	b.x     =  a.yz;
	b.y     =  a.zx;
	b.z     =  a.xy;

	b.xy   =   a.z;
	b.zx   =   a.y;
	b.yz   =   a.x;

	b.xyz  =  a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_UnderBar(GA3Ezx a)   // in 3D, OverBar and UnderBar coincide
// Table 4.3 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA3Ezx b;		// b blade wedge a blade = pseudovector blade
// UnderBar(Blade[i]) = xyz*Reciprocal(Blade[i]) ;

	b.q     =  a.xyz;	// xyz

	b.x     =  a.yz;
	b.y     =  a.zx;
	b.z     =  a.xy;

	b.xy   =   a.z;
	b.zx   =   a.y;
	b.yz   =   a.x;

	b.xyz  =  a.q;

	return b;
}


//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Reverse(GA3Ezx w)	
// (A,  B, C, D,  -E,-F,-G, -H) 
{
	GA3Ezx v;
	v.q =  w.q;

	v.x = w.x;
	v.y = w.y;
	v.z = w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Hermitian(GA3Ezx w)	
//(A,  B, C, D, -E,-F,-G, -H)  
// correct for both 2x2 complex matrix and real 4x4 matrix implementations
{
	GA3Ezx v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Transpose(GA3Ezx w)	
//(A,  B, C, D, -E,-F,-G, -H)  
{
	GA3Ezx v;

	v.q =  w.q;

	v.x =  w.x;
	v.y =  w.y;
	v.z =  w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Parity(GA3Ezx w)	
//(A,  B, C, D, -E,-F,-G, -H)  
{
	GA3Ezx v;

	v.q =  w.q;

	v.x =  -w.x;
	v.y =  -w.y;
	v.z =  -w.z;

	v.xy = w.xy;
	v.zx = w.zx;
	v.yz = w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Conjugation(GA3Ezx w)
// (A, -B,-C,-D, -E,-F,-G, -H)
{
	GA3Ezx v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = -w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_CliffordConjugation(GA3Ezx w)
// (A, -B,-C,-D, -E,-F,-G,  H) 
{
	GA3Ezx v;
	v.q =  w.q;

	v.x = -w.x;
	v.y = -w.y;
	v.z = -w.z;

	v.xy = -w.xy;
	v.zx = -w.zx;
	v.yz = -w.yz;

	v.xyz = w.xyz;

	return v;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Dual(GA3Ezx w)   // return w*I_inv 
// Dual(r) = u = (h, g,-f, e, -d, c,-b, -a)

{
	GA3Ezx v;
//	GA3Ezx I_inv;
//	I_inv.xyz = -1;
//	v = GA3Ezx_Product(w,I_inv);

	v.q =  w.xyz;

	v.x =  w.yz;
	v.y = -w.zx;
	v.z =  w.xy;

	v.xy = -w.z;
	v.zx =  w.y;
	v.yz = -w.x;

	v.xyz = -w.q;

	return v;
}



//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_DorstDual(GA3Ezx a)  //DorstDual(r) = u = (h, g,f,e, -d,-c,-b, -a)
{
	GA3Ezx b;
//	GA3Ezx I_inv;	I_inv.xyz = -1;

//	b = GA3Ezx_LeftContraction(a,I_inv);

//	b = LeftContraction(a,I_inv); 

	b.q =  a.xyz;

	b.x =  a.yz;
	b.y =  a.zx;
	b.z =  a.xy;

	b.xy = -a.z;
	b.zx = -a.y;
	b.yz = -a.x;

	b.xyz = -a.q;

	return b;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_DorstUnDual(GA3Ezx a)  //DorstUnDual(r) = u = (-h,-g,-f,-e,  d, c, b,  a)
{
	GA3Ezx b;
//	GA3Ezx I;	I.xyz = 1;

//	b = GA3Ezx_LeftContraction(a,I);

	b.q = -a.xyz;

	b.x = -a.yz;
	b.y = -a.zx;
	b.z = -a.xy;

	b.xy =  a.z;
	b.zx =  a.y;
	b.yz =  a.x;

	b.xyz =  a.q;

	return b;
}


//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Add(GA3Ezx u,GA3Ezx v)
{
	GA3Ezx w;
	w.q =   u.q   + v.q  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.z   = u.z   + v.z  ;
	w.xy  = u.xy  + v.xy ;
	w.zx  = u.zx  + v.zx ;
	w.yz  = u.yz  + v.yz ;
	w.xyz = u.xyz + v.xyz;

	return w;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Subtract(GA3Ezx u, GA3Ezx v)
{
	GA3Ezx w;
	w.q =   u.q   - v.q  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.z   = u.z   - v.z  ;
	w.xy  = u.xy  - v.xy ;
	w.zx  = u.zx  - v.zx ;
	w.yz  = u.yz  - v.yz ;
	w.xyz = u.xyz - v.xyz;

	return w;
}


//////////////////////////////////////////////////////

int GA3Ezx_Equal(GA3Ezx u, GA3Ezx v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.x ==v.x )&&
			(u.y ==v.y )&&
			(u.z ==v.z )&&

			(u.xy==v.xy)&&
			(u.zx==v.zx)&&
			(u.yz==v.yz)&&

			(u.xyz==v.xyz);
	return result;
}

//////////////////////////////////////////////////////

int GA3Ezx_Not_Equal(GA3Ezx u, GA3Ezx v)
{
	int result;
	result = 	(u.q !=v.q )||

			(u.x !=v.x )||
			(u.y !=v.y )||
			(u.z !=v.z )||

			(u.xy!=v.xy)||
			(u.zx!=v.zx)||
			(u.yz!=v.yz)||

			(u.xyz!=v.xyz);
	return result;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Product(GA3Ezx a, GA3Ezx b)
{

/*

 *    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     q     xy   -zx    y    -z     xyz   yz   
 y    |  y    -xy    q     yz   -x     xyz   z     zx   
 z    |  z     zx   -yz    q     xyz   x    -y     xy   
 xy   |  xy   -y     x     xyz  -q     yz   -zx   -z    
 zx   |  zx    z     xyz  -x    -yz   -q     xy   -y    
 yz   |  yz    xyz  -z     y     zx   -xy   -q    -x    
 xyz  |  xyz   yz    zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

	c.q   = a.q*b.q   + a.x*b.x   + a.y*b.y   + a.z*b.z   - a.xy*b.xy  - a.zx*b.zx  - a.yz*b.yz  - a.xyz*b.xyz;
	c.x   = a.q*b.x   + a.x*b.q   - a.y*b.xy  + a.z*b.zx  + a.xy*b.y   - a.zx*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y   = a.q*b.y   + a.x*b.xy  + a.y*b.q   - a.z*b.yz  - a.xy*b.x   - a.zx*b.xyz + a.yz*b.z   - a.xyz*b.zx ;
	c.z   = a.q*b.z   - a.x*b.zx  + a.y*b.yz  + a.z*b.q   - a.xy*b.xyz + a.zx*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xy  = a.q*b.xy  + a.x*b.y   - a.y*b.x   + a.z*b.xyz + a.xy*b.q   + a.zx*b.yz  - a.yz*b.zx  + a.xyz*b.z  ;
	c.zx  = a.q*b.zx  - a.x*b.z   + a.y*b.xyz + a.z*b.x   - a.xy*b.yz  + a.zx*b.q   + a.yz*b.xy  + a.xyz*b.y  ;
	c.yz  = a.q*b.yz  + a.x*b.xyz + a.y*b.z   - a.z*b.y   + a.xy*b.zx  - a.zx*b.xy  + a.yz*b.q   + a.xyz*b.x  ;
	c.xyz = a.q*b.xyz + a.x*b.yz  + a.y*b.zx  + a.z*b.xy  + a.xy*b.z   + a.zx*b.y   + a.yz*b.x   + a.xyz*b.q  ;

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Multiply_By_Constant(GA3Ezx u, double a)
{
	GA3Ezx w;

	w.q =   u.q*a;

	w.x   = u.x*a;
	w.y   = u.y*a;
	w.z   = u.z*a;

	w.xy  = u.xy*a;
	w.zx  = u.zx*a;
	w.yz  = u.yz*a;

	w.xyz = u.xyz*a;

	return w;
}


//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Divide_By_Constant(GA3Ezx u, double a)
{
	GA3Ezx w;

	w.q =   u.q/a;

	w.x   = u.x/a;
	w.y   = u.y/a;
	w.z   = u.z/a;

	w.xy  = u.xy/a;
	w.zx  = u.zx/a;
	w.yz  = u.yz/a;

	w.xyz = u.xyz/a;

	return w;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Wedge(GA3Ezx a, GA3Ezx b) {

	GA3Ezx c;

/*

 ^    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     0     xy   -zx    0     0     xyz   0    
 y    |  y    -xy    0     yz    0     xyz   0     0    
 z    |  z     zx   -yz    0     xyz   0     0     0    
 xy   |  xy    0     0     xyz   0     0     0     0    
 zx   |  zx    0     xyz   0     0     0     0     0    
 yz   |  yz    xyz   0     0     0     0     0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/

	c.q   =  a.q*b.q;  
	c.x   =  a.q*b.x   + a.x*b.q;  
	c.y   =  a.q*b.y   + a.y*b.q;  
	c.z   =  a.q*b.z   + a.z*b.q; 
	c.xy  =  a.q*b.xy  + a.x*b.y  - a.y*b.x  + a.xy*b.q;  
	c.zx  =  a.q*b.zx  - a.x*b.z  + a.z*b.x  + a.zx*b.q;  
	c.yz  =  a.q*b.yz  + a.y*b.z  - a.z*b.y  + a.yz*b.q;  
	c.xyz =  a.q*b.xyz + a.x*b.yz + a.y*b.zx + a.xy*b.z + a.z*b.xy + a.zx*b.y + a.yz*b.x + a.xyz*b.q;

	return c;

}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_AntiWedge(GA3Ezx a, GA3Ezx b)
{
	GA3Ezx c;
/*

Lengyel's AntiWedge
 V    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0     q     0     y    
 z    |  0     0     0     0     q     0     0     z    
 xy   |  0     0     0     q     0    -x     y     xy   
 zx   |  0     0     q     0     x     0    -z     zx   
 yz   |  0     q     0     0    -y     z     0     yz   
 xyz  |  q     x     y     z     xy    zx    yz    xyz  

*/

//	c = OverBar(Wedge(UnderBar(a),UnderBar(b)));

// AntiWedge
	c.q   =  + a.q  *b.xyz + a.x  *b.yz  + a.y  *b.zx  + a.z  *b.xy
                 + a.xy *b.z   + a.zx *b.y   + a.yz *b.x   + a.xyz*b.q   ; 
	c.x   =  + a.x  *b.xyz - a.xy *b.zx  + a.zx *b.xy  + a.xyz*b.x   ; 
	c.y   =  + a.y  *b.xyz + a.xy *b.yz  - a.yz *b.xy  + a.xyz*b.y   ; 
	c.z   =  + a.z  *b.xyz - a.zx *b.yz  + a.yz *b.zx  + a.xyz*b.z   ; 
	c.xy  =  + a.xy *b.xyz + a.xyz*b.xy  ; 
	c.zx  =  + a.zx *b.xyz + a.xyz*b.zx  ; 
	c.yz  =  + a.yz *b.xyz + a.xyz*b.yz  ; 
	c.xyz =  + a.xyz*b.xyz ; 

	return c;
}

////////////////////////////////////////////////////// 


GA3Ezx GA3Ezx_Regressive(GA3Ezx a, GA3Ezx b) 
{
	GA3Ezx c;	

/*
Hestenes' Regressive (differs from AntiWedge in signs)

 V    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     q    
 x    |  0     0     0     0     0     0     q     x    
 y    |  0     0     0     0     0     q     0     y    
 z    |  0     0     0     0     q     0     0     z    
 xy   |  0     0     0     q     0     x    -y     xy   
 zx   |  0     0     q     0    -x     0     z     zx   
 yz   |  0     q     0     0     y    -z     0     yz   
 xyz  |  q     x     y     z     xy    zx    yz    xyz  

*/

//	GA3Ezx I, I_inv;
//	I.xyz = 1;
//	I_inv.xyz = -1;

//	c = ((a*I_inv)^(b*I_inv))*I;

	c.q = a.q*b.xyz + a.x*b.yz + a.y*b.zx + a.z*b.xy + a.xy*b.z + a.zx*b.y + a.yz*b.x + a.xyz*b.q ;

	c.x = + ( a.x*b.xyz + a.xy*b.zx - a.zx*b.xy + a.xyz*b.x) ;
	c.y = + ( a.y*b.xyz - a.xy*b.yz + a.yz*b.xy + a.xyz*b.y) ;
	c.z = + ( a.z*b.xyz + a.zx*b.yz - a.yz*b.zx + a.xyz*b.z) ;

	c.xy = + ( a.xy*b.xyz + a.xyz*b.xy) ;
	c.zx = + ( a.zx*b.xyz + a.xyz*b.zx) ;
	c.yz = + ( a.yz*b.xyz + a.xyz*b.yz) ;

	c.xyz = + a.xyz*b.xyz ;

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_RegressiveViaFormula(GA3Ezx a, GA3Ezx b)
{
	GA3Ezx c;

	GA3Ezx Ixyz,I_inv;
	
	Ixyz = GA3Ezx_Zero();	I_inv = GA3Ezx_Zero();
	Ixyz.xyz = 1;	I_inv.xyz = -1;

	c = GA3Ezx_Product(GA3Ezx_Wedge(GA3Ezx_Product(a,I_inv),GA3Ezx_Product(b,I_inv)),Ixyz);

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_LowerRightViaFormula(GA3Ezx a, GA3Ezx b)  // Like a mirror of Wedge along rising diagonal
{

/*


LowerRightViaFormula(Blade[i],Blade[j])
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     xyz  
 x    |  0     0     0     0     0     0     xyz   yz   
 y    |  0     0     0     0     0    -xyz   0    -zx   
 z    |  0     0     0     0     xyz   0     0     xy   
 xy   |  0     0     0     xyz   0    -yz    zx   -z    
 zx   |  0     0    -xyz   0     yz    0    -xy    y    
 yz   |  0     xyz   0     0    -zx    xy    0    -x    
 xyz  |  xyz   yz   -zx    xy   -z     y    -x    -q    


*/
// Wedge(a*I_inv,I*b);

	GA3Ezx c;

	GA3Ezx Ixyz,I_inv;
	
	Ixyz = GA3Ezx_Zero();	I_inv = GA3Ezx_Zero();
	Ixyz.xyz = 1;	I_inv.xyz = -1;

	c = GA3Ezx_Wedge(GA3Ezx_Product(a,I_inv),GA3Ezx_Product(Ixyz,b));
	

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Expander(GA3Ezx a, GA3Ezx b)
{ 
/*
Terms with increased rank
 >    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     0     xy   -zx    0     0     xyz   0    
 y    |  0    -xy    0     yz    0     xyz   0     0    
 z    |  0     zx   -yz    0     xyz   0     0     0    
 xy   |  0     0     0     xyz   0     0     0     0    
 zx   |  0     0     xyz   0     0     0     0     0    
 yz   |  0     xyz   0     0     0     0     0     0    
 xyz  |  0     0     0     0     0     0     0     0    

*/
	GA3Ezx c;

	c.q    =  0 ; 

	c.x    =  0 ; 
	c.y    =  0 ; 
	c.z    =  0 ; 

	c.xy   =  + a.x*b.y - a.y*b.x    ; 
	c.zx   =  - a.x*b.z + a.z*b.x    ; 
	c.yz   =  + a.y*b.z - a.z*b.y    ; 

	c.xyz  =  + a.x*b.yz + a.y*b.zx + a.z*b.xy + a.xy*b.z + a.zx*b.y + a.yz*b.x    ; 

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Conserver(GA3Ezx a, GA3Ezx b)
{  

/*

Terms with preserved rank
 =    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     0     0     0     0     0     0     0    
 y    |  y     0     0     0     0     0     0     0    
 z    |  z     0     0     0     0     0     0     0    
 xy   |  xy    0     0     0     0     yz   -zx    0    
 zx   |  zx    0     0     0    -yz    0     xy    0    
 yz   |  yz    0     0     0     zx   -xy    0     0    
 xyz  |  xyz   0     0     0     0     0     0     0    

*/
	GA3Ezx c;

	c.q    =  + a.q*b.q    ; 

	c.x    =  + a.q*b.x    + a.x*b.q    ; 
	c.y    =  + a.q*b.y    + a.y*b.q    ; 	
	c.z    =  + a.q*b.z    + a.z*b.q    ; 

	c.xy   =  + a.q*b.xy   + a.xy*b.q    + a.zx*b.yz   - a.yz*b.zx ; 
	c.zx   =  + a.q*b.zx   - a.xy*b.yz   + a.zx*b.q    + a.yz*b.xy ; 
	c.yz   =  + a.q*b.yz   + a.xy*b.zx   - a.zx*b.xy   + a.yz*b.q  ; 

	c.xyz  =  + a.q*b.xyz  + a.xyz *b.q  ; 

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Shrinker(GA3Ezx a, GA3Ezx b)
{

/*

Terms with reduced rank
 <    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     y    -z     0     yz   
 y    |  0     0     q     0    -x     0     z     zx   
 z    |  0     0     0     q     0     x    -y     xy   
 xy   |  0    -y     x     0    -q     0     0    -z    
 zx   |  0     z     0    -x     0    -q     0    -y    
 yz   |  0     0    -z     y     0     0    -q    -x    
 xyz  |  0     yz    zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

// Shrinker equation set
	c.q    =  + a.x*b.x + a.y*b.y + a.z*b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz *b.xyz; 

	c.x    =  - a.y*b.xy + a.z*b.zx + a.xy*b.y - a.zx*b.z - a.yz*b.xyz - a.xyz*b.yz ; 
	c.y    =  + a.x*b.xy - a.z*b.yz - a.xy*b.x - a.zx*b.xyz + a.yz*b.z - a.xyz*b.zx ; 
	c.z    =  - a.x*b.zx + a.y*b.yz - a.xy*b.xyz + a.zx*b.x - a.yz*b.y - a.xyz*b.xy ; 

	c.xy   =  + a.z*b.xyz + a.xyz*b.z  ; 
	c.zx   =  + a.y*b.xyz + a.xyz*b.y ; 
	c.yz   =  + a.x*b.xyz + a.xyz*b.x ; 

	c.xyz  =  0  ; 

	return c;
}

///////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Symmetric(GA3Ezx a, GA3Ezx b)
{ 
/*

Symmetric Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  x     q     0     0     0     0     xyz   yz   
 y    |  y     0     q     0     0     xyz   0     zx   
 z    |  z     0     0     q     xyz   0     0     xy   
 xy   |  xy    0     0     xyz  -q     0     0    -z    
 zx   |  zx    0     xyz   0     0    -q     0    -y    
 yz   |  yz    xyz   0     0     0     0    -q    -x    
 xyz  |  xyz   yz    zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

//	c = (a*b + b*a)/2;

	c.q   =  + a.q*b.q   + a.x*b.x   + a.y *b.y   + a.z  *b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz*b.xyz ; 
	c.x   =  + a.q*b.x   + a.x*b.q   - a.yz*b.xyz - a.xyz*b.yz  ; 
	c.y   =  + a.q*b.y   + a.y*b.q   - a.zx*b.xyz - a.xyz*b.zx  ; 
	c.z   =  + a.q*b.z   + a.z*b.q   - a.xy*b.xyz - a.xyz*b.xy  ; 
	c.xy  =  + a.q*b.xy  + a.z*b.xyz + a.xy*b.q   + a.xyz*b.z   ; 
	c.zx  =  + a.q*b.zx  + a.y*b.xyz + a.zx*b.q   + a.xyz*b.y   ; 
	c.yz  =  + a.q*b.yz  + a.x*b.xyz + a.yz*b.q   + a.xyz*b.x   ; 
	c.xyz =  + a.q*b.xyz + a.x*b.yz  + a.y *b.zx  + a.z*b.xy  + a.xy*b.z + a.zx*b.y + a.yz*b.x + a.xyz*b.q   ; 

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_AntiSymmetric(GA3Ezx a, GA3Ezx b)
{ 
/*

AntiSymmetric Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     0     xy   -zx    y    -z     0     0    
 y    |  0    -xy    0     yz   -x     0     z     0    
 z    |  0     zx   -yz    0     0     x    -y     0    
 xy   |  0    -y     x     0     0     yz   -zx    0    
 zx   |  0     z     0    -x    -yz    0     xy    0    
 yz   |  0     0    -z     y     zx   -xy    0     0    
 xyz  |  0     0     0     0     0     0     0     0    

*/
	GA3Ezx c;

//	c = (a*b - b*a)/2;

	c.q   =  0; 
	c.x   =  - a.y*b.xy  + a.z*b.zx  + a.xy*b.y   - a.zx*b.z   ; 
	c.y   =  + a.x*b.xy  - a.z*b.yz  - a.xy*b.x   + a.yz*b.z   ; 
	c.z   =  - a.x*b.zx  + a.y*b.yz  + a.zx*b.x   - a.yz*b.y   ; 
	c.xy  =  + a.x*b.y   - a.y*b.x   + a.zx*b.yz  - a.yz*b.zx  ; 
	c.zx  =  - a.x*b.z   + a.z*b.x   - a.xy*b.yz  + a.yz*b.xy  ; 
	c.yz  =  + a.y*b.z   - a.z*b.y   + a.xy*b.zx  - a.zx*b.xy  ; 
	c.xyz =  0; 

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Inner(GA3Ezx a, GA3Ezx b) // Based Upon Bromborsky (A|B)
{ 
/*

Inner Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  0     0     0     0     0     0     0     0    
 x    |  0     q     0     0     y    -z     0     yz   
 y    |  0     0     q     0    -x     0     z     zx   
 z    |  0     0     0     q     0     x    -y     xy   
 xy   |  0    -y     x     0    -q     0     0    -z    
 zx   |  0     z     0    -x     0    -q     0    -y    
 yz   |  0     0    -z     y     0     0    -q    -x    
 xyz  |  0     yz   -zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

	c.q =  + a.x*b.x  + a.y*b.y  + a.z*b.z  - a.xy*b.xy  - a.zx*b.zx  - a.yz*b.yz  - a.xyz*b.xyz ;
	c.x =             - a.y*b.xy + a.z*b.zx + a.xy*b.y   - a.zx*b.z   - a.yz*b.xyz - a.xyz*b.yz ;
	c.y =  + a.x*b.xy            - a.z*b.yz - a.xy*b.x   - a.zx*b.xyz + a.yz*b.z   - a.xyz*b.zx ;
	c.z =  - a.x*b.zx + a.y*b.yz            - a.xy*b.xyz + a.zx*b.x   - a.yz*b.y   - a.xyz*b.xy ;
	c.xy = + a.z*b.xyz + a.xyz*b.z  ;
	c.zx = + a.y*b.xyz + a.xyz*b.y  ;
	c.yz = + a.x*b.xyz + a.xyz*b.x ;
	c.xyz = 0;

	return c;
}

//////////////////////////////////////////////////////  // Based Upon Bromborsky (A<B)

GA3Ezx GA3Ezx_LeftContraction (GA3Ezx a, GA3Ezx b)
{ 
/*

LeftContraction Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     x     y     z     xy    zx    yz    xyz  
 x    |  0     q     0     0     y    -z     0     yz   
 y    |  0     0     q     0    -x     0     z     zx   
 z    |  0     0     0     q     0     x    -y     xy   
 xy   |  0     0     0     0    -q     0     0    -z    
 zx   |  0     0     0     0     0    -q     0    -y    
 yz   |  0     0     0     0     0     0    -q    -x    
 xyz  |  0     0     0     0     0     0     0    -q    

*/
	GA3Ezx c;

	c.q = a.q*b.q + a.x*b.x + a.y*b.y + a.z*b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz*b.xyz ;

	c.x = + (a.q*b.x - a.y*b.xy + a.z*b.zx - a.yz*b.xyz) ;
	c.y = + (a.q*b.y + a.x*b.xy - a.z*b.yz - a.zx*b.xyz) ;
	c.z = + (a.q*b.z - a.x*b.zx + a.y*b.yz - a.xy*b.xyz) ;

	c.xy = + (a.q*b.xy + a.z*b.xyz) ;
	c.zx = + (a.q*b.zx + a.y*b.xyz) ;
	c.yz = + (a.q*b.yz + a.x*b.xyz) ;

	c.xyz = + a.q*b.xyz ;

	return c;
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_RightContraction (GA3Ezx a, GA3Ezx b)  // Based Upon Bromborsky (A>B)
{ 
/*

RightContraction Product
 ?    |  q     x     y     z     xy    zx    yz    xyz  
-------------------------------------------------------
 q    |  q     0     0     0     0     0     0     0    
 x    |  x     q     0     0     0     0     0     0    
 y    |  y     0     q     0     0     0     0     0    
 z    |  z     0     0     q     0     0     0     0    
 xy   |  xy   -y     x     0    -q     0     0     0    
 zx   |  zx    z     0    -x     0    -q     0     0    
 yz   |  yz    0    -z     y     0     0    -q     0    
 xyz  |  xyz   yz   -zx    xy   -z    -y    -x    -q    

*/
	GA3Ezx c;

	c.q = a.q*b.q + a.x*b.x + a.y*b.y + a.z*b.z - a.xy*b.xy - a.zx*b.zx - a.yz*b.yz - a.xyz*b.xyz ;

	c.x = + ( a.x*b.q + a.xy*b.y - a.zx*b.z - a.xyz*b.yz) ;
	c.y = + ( a.y*b.q - a.xy*b.x + a.yz*b.z - a.xyz*b.zx) ;
	c.z = + ( a.z*b.q + a.zx*b.x - a.yz*b.y - a.xyz*b.xy) ;

	c.xy =  ( a.xy*b.q + a.xyz*b.z) ;
	c.zx =  ( a.zx*b.q - a.xyz*b.y) ; 
	c.yz =  ( a.yz*b.q + a.xyz*b.x) ;

	c.xyz = + a.xyz*b.q ;

	return c;
}

//////////////////////////////////////////////////////

double GA3Ezx_Determinant(GA3Ezx A) {

	GA3Ezx B, C;
	double a, b, det;


	B = GA3Ezx_CliffordConjugation(A);
	C = GA3Ezx_Product(B,A);

	a = C.q;
	b = C.xyz;

	det = (a*a + b*b);

	return(det);
}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Adjugate(GA3Ezx a)
{

	GA3Ezx u;

//Product(CliffordConjugation(a),Conjugation(Product(a,CliffordConjugation(a)))) = Adjugate(a);

	u = GA3Ezx_Product(GA3Ezx_CliffordConjugation(a),GA3Ezx_Conjugation(GA3Ezx_Product(a,GA3Ezx_CliffordConjugation(a))));

	return u;

}

//////////////////////////////////////////////////////

GA3Ezx GA3Ezx_Reciprocal(GA3Ezx a)
{
	double b;
	GA3Ezx c,d;

	b = GA3Ezx_Determinant(a);

	c = GA3Ezx_Adjugate(a);

	d.q = (c.q/b);	// divide by determinant

	d.x = (c.x/b);
	d.y = (c.y/b);
	d.z = (c.z/b);

	d.xy = (c.xy/b);
	d.zx = (c.zx/b);
	d.yz = (c.yz/b);

	d.xyz = (c.xyz/b);

	return d;
}

////////////////////////////////////////////////////// 

GA3Ezx	GA3Ezx_Magic(GA3Ezx u, int i) {

	GA3Ezx MV;
	double a, b,c,d, e,f,g, h;

	a = u.q;	b = u.x;	c = u.y;	d = u.z;
	e = u.xy;	f = u.zx;	g = u.yz;	h = u.xyz;

	if(i == 0) MV = GA3Ezx_Set( a, +b,+c,+d, +e,+f,+g, +h) ;
	if(i == 1) MV = GA3Ezx_Set( a, +b,+d,+c, -f,-e,-g, -h) ;
	if(i == 2) MV = GA3Ezx_Set( a, +c,+b,+d, -e,-g,-f, -h) ; 
	if(i == 3) MV = GA3Ezx_Set( a, +c,+d,+b, +g,+e,+f, +h) ; 
	if(i == 4) MV = GA3Ezx_Set( a, +d,+b,+c, +f,+g,+e, +h) ; 
	if(i == 5) MV = GA3Ezx_Set( a, +d,+c,+b, -g,-f,-e, -h) ;
	
	return MV;
}

//////////////////////////////////////////////////////

GA3Ezx	GA3Ezx_Comp(GA3Ezx u, int i) {

	GA3Ezx MV;
	double a, b,c,d, e,f,g, h;
	int Mask, Sign[8];

	Mask = 16;
	if((Mask & i) == 0) Sign[4] = +1; else Sign[4] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[3] = +1; else Sign[3] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[2] = +1; else Sign[2] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[1] = +1; else Sign[1] = -1;  Mask >>= 1;
	if((Mask & i) == 0) Sign[0] = +1; else Sign[0] = -1;  Mask >>= 1;
	Sign[5] = Sign[2]*Sign[1]*Sign[0];
	Sign[6] = Sign[3]*Sign[1]*Sign[0];
	Sign[7] = Sign[4]*Sign[1]*Sign[0];

	a = Sign[4]*u.q;   b = Sign[3]*u.x;   c = Sign[2]*u.y;   d = Sign[1]*u.z;
	e = Sign[0]*u.xy;  f = Sign[5]*u.zx;  g = Sign[6]*u.yz;  h = Sign[7]*u.xyz;

	MV = GA3Ezx_Set(a, b,c,d, e,f,g, h) ; 
	
	return MV;

}


//////////////////////////////////////////////////////

//	int main(void) { return 0; }
