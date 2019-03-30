// Routines for Geometric Algebra in Five Dimensional metric(4,1) space
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

//////////////////////////////////////////////////////

typedef struct
{
	long q ;
	long w,x,y,z,t ;
	long wx,wy,wz,wt,xy,xz,xt,yz,yt,zt ;
	long wxy,wxz,wxt,wyz,wyt,wzt,xyz,xyt,xzt,yzt ;
	long wxyz,wxyt,wxzt,wyzt,xyzt; 
	long wxyzt;
} GA5_4_1;


//////////////////////////////////////////////////////

typedef struct
{
	long r,i ;
} complex;


//////////////////////////////////////////////////////

GA5_4_1 Zero(void) 
{

	GA5_4_1 a;
	
 	a.q = 0.0 ;
	a.w = 0.0 ; a.x = 0.0 ; a.y = 0.0 ; a.z = 0.0 ; a.t = 0.0 ; 
	a.wx = 0.0 ; a.wy = 0.0 ; a.wz = 0.0 ; a.wt = 0.0 ; a.xy = 0.0 ;
	a.xz = 0.0 ; a.xt = 0.0 ; a.yz = 0.0 ; a.yt = 0.0 ; a.zt = 0.0 ; 
	a.wxy = 0.0 ; a.wxz = 0.0 ; a.wxt = 0.0 ; a.wyz = 0.0 ; a.wyt = 0.0 ;
	a.wzt = 0.0 ; a.xyz = 0.0 ; a.xyt = 0.0 ; a.xzt = 0.0 ; a.yzt = 0.0 ;   
	a.wxyz = 0.0 ; a.wxyt = 0.0 ; a.wxzt = 0.0 ; a.wyzt = 0.0 ; a.xyzt = 0.0 ;   
	a.wxyzt = 0.0 ; 

	return a;

}


//////////////////////////////////////////////////////


GA5_4_1 Set(long q, 
		long w, long x, long y, long z, long t, 
		long wx, long wy, long wz, long wt, long xy, 
		long xz, long xt, long yz, long yt, long zt, 
		long wxy, long wxz, long wxt, long wyz, long wyt, 
		long wzt, long xyz, long xyt, long xzt, long yzt, 
		long wxyz, long wxyt, long wxzt, long wyzt, long xyzt, 
		long wxyzt) 
{

	GA5_4_1 a;
	
 	a.q = q ;
	a.w = w ; a.x = x ; a.y = y ; a.z = z ; a.t = t ; 
	a.wx = wx ; a.wy = wy ; a.wz = wz ; a.wt = wt ; a.xy = xy ;
	a.xz = xz ; a.xt = xt ; a.yz = yz ; a.yt = yt ; a.zt = zt ; 
	a.wxy = wxy ; a.wxz = wxz ; a.wxt = wxt ; a.wyz = wyz ; a.wyt = wyt ;
	a.wzt = wzt ; a.xyz = xyz ; a.xyt = xyt ; a.xzt = xzt ; a.yzt = yzt ;   
	a.wxyz = wxyz ; a.wxyt = wxyt ; a.wxzt = wxzt ; a.wyzt = wyzt ; a.xyzt = xyzt ;   
	a.wxyzt = wxyzt ; 

	return a;
}

//////////////////////////////////////////////////////

// Necessary forward declarations

GA5_4_1 LeftContraction ( GA5_4_1 a,  GA5_4_1 b) ;
GA5_4_1 Product( GA5_4_1 a,  GA5_4_1 b) ;

//////////////////////////////////////////////////////

void PrintlnMV(GA5_4_1 v) {

	printf("(%10ld, \n", v.q) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld, \n", v.w, v.x, v.y, v.z, v.t) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld, \n", 
		v.wx, v.wy, v.wz, v.wt, v.xy, v.xz, v.xt, v.yz, v.yt, v.zt) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld, \n", 
		v.wxy, v.wxz, v.wxt, v.wyz, v.wyt, v.wzt, v.xyz, v.xyt, v.xzt, v.yzt) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld, \n", v.wxyz, v.wxyt, v.wxzt, v.wyzt, v.xyzt) ;
	printf("%10ld) \n", v.wxyzt);

}


//////////////////////////////////////////////////////

void PrintMV(GA5_4_1 v) {

	printf("(%10ld, ", v.q) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld, ", v.w, v.x, v.y, v.z, v.t) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld, ", 
		v.wx, v.wy, v.wz, v.wt, v.xy, v.xz, v.xt, v.yz, v.yt, v.zt) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld,%10ld, ", 
		v.wxy, v.wxz, v.wxt, v.wyz, v.wyt, v.wzt, v.xyz, v.xyt, v.xzt, v.yzt) ;
	printf("%10ld,%10ld,%10ld,%10ld,%10ld, ", v.wxyz, v.wxyt, v.wxzt, v.wyzt, v.xyzt) ;
	printf("%10ld) \n", v.wxyzt);

}

//////////////////////////////////////////////////////

GA5_4_1 OverBar(GA5_4_1 a)
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA5_4_1 b;		// a blade wedge b blade = pseudovector blade wxyzt
// a^OverBar(a) = pseudovector blade wxyzt

	b.q = a.wxyzt;

	b.w =  a.xyzt; 	b.x = -a.wyzt; 	b.y =  a.wxzt; 	b.z = -a.wxyt; 	b.t =  a.wxyz;

	b.wx =  a.yzt; 	b.wy = -a.xzt; 	b.wz =  a.xyt; 	b.wt = -a.xyz; 	b.xy =  a.wzt;
	b.xz = -a.wyt; 	b.xt =  a.wyz; 	b.yz =  a.wxt; 	b.yt = -a.wxz; 	b.zt =  a.wxy;

	b.wxy =  a.zt; 	b.wxz = -a.yt; 	b.wxt =  a.yz; 	b.wyz =  a.xt; 	b.wyt = -a.xz;
	b.wzt =  a.xy; 	b.xyz = -a.wt; 	b.xyt =  a.wz; 	b.xzt = -a.wy; 	b.yzt =  a.wx; 

	b.wxyz =  a.t; 	b.wxyt = -a.z; 	b.wxzt =  a.y; 	b.wyzt = -a.x; 	b.xyzt =  a.w;

	b.wxyzt = a.q;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 UnderBar(GA5_4_1 a)   
// Table 4.4 Lengyel, Volume 1 Mathematics, Foundations of Game Engine Development
{
	GA5_4_1 b;		// b blade wedge a blade = pseudovector blade


	b.q = a.wxyzt;

	b.w =  a.xyzt; 	b.x = -a.wyzt; 	b.y =  a.wxzt; 	b.z = -a.wxyt; 	b.t =  a.wxyz;

	b.wx =  a.yzt; 	b.wy = -a.xzt; 	b.wz =  a.xyt; 	b.wt = -a.xyz; 	b.xy =  a.wzt;
	b.xz = -a.wyt; 	b.xt =  a.wyz; 	b.yz =  a.wxt; 	b.yt = -a.wxz; 	b.zt =  a.wxy;

	b.wxy =  a.zt; 	b.wxz = -a.yt; 	b.wxt =  a.yz; 	b.wyz =  a.xt; 	b.wyt = -a.xz;
	b.wzt =  a.xy; 	b.xyz = -a.wt; 	b.xyt =  a.wz; 	b.xzt = -a.wy; 	b.yzt =  a.wx; 

	b.wxyz =  a.t; 	b.wxyt = -a.z; 	b.wxzt =  a.y; 	b.wyzt = -a.x; 	b.xyzt =  a.w;

	b.wxyzt = a.q;

	return b;
}  


//////////////////////////////////////////////////////

GA5_4_1 Reverse(GA5_4_1 a)	// change sign based upon multivector product reversal.
// Ex:  xy -> -yx, xyz -> -zyx, xyzt -> tzyx, wxyzt => tzyxw
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w =  a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = a.t;
	b.wx = - a.wx;		b.wy = - a.wy;		b.wz = - a.wz;		b.wt = - a.wt;		b.xy = - a.xy;
	b.xz = - a.xz;		b.yz = - a.yz;		b.xt = - a.xt;		b.yt = - a.yt;		b.zt = - a.zt;
	b.wxy = - a.wxy;	b.wxz = - a.wxz;	b.wxt = - a.wxt;	b.wyz = - a.wyz;	b.wyt = - a.wyt;
	b.wzt = - a.wzt;	b.xyz = - a.xyz;	b.xyt = - a.xyt;	b.xzt = - a.xzt;	b.yzt = - a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt =  a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Reverse_Vector_Quad(GA5_4_1 a)	// change sign based upon multivector product reversal.
// Ex:  xy -> -yx, xyz -> -zyx, xyzt -> tzyx, wxyzt => tzyxw
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = a.wxy;		b.wxz = a.wxz;		b.wxt = a.wxt;		b.wyz = a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt =  a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Involution(GA5_4_1 a)	
// Corresponds to WPT parity transform: w -> -w,  x -> -x, y -> -y, z ->-z, t -> -t
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = -a.wxt;		b.wyz = -a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = -a.xyt;		b.xzt = -a.xzt;		b.yzt = -a.yzt;
	b.wxyz =   a.wxyz;	b.wxyt =   a.wxyt;	b.wxzt =   a.wxzt;	b.wyzt =   a.wyzt;	b.xyzt =   a.xyzt;
	b.wxyzt = -a.wxyzt;

	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 Transpose(GA5_4_1 a)	// verified
//	+ a*q
//	- b*w + c*x + d*y + e*z - f*t
//	+ g*wx + h*wy + j*wz - k*wt - l*xy - m*xz + n*xt - p*yz + r*yt + s*zt
//	+ S*wxy + R*wxz - P*wxt + N*wyz - M*wyt - L*wzt - K*xyz + J*xyt + H*xzt + G*yzt
//	- F*wxyz + E*wxyt + D*wxzt + C*wyzt - B*xyzt
//	+ A*wxyzt 
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = - a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = - a.wt;		b.xy = - a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = a.wxy;		b.wxz = a.wxz;		b.wxt = -a.wxt;		b.wyz = a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = a.wxyzt;

	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 Conjugation(GA5_4_1 a)	
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t =  -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt =  -a.wt;		b.xy =  -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = -a.xt;		b.yt = -a.yt;		b.zt = -a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = -a.wxt;		b.wyz = -a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = -a.xyt;		b.xzt = -a.xzt;		b.yzt = -a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = -a.wxyzt;

	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 CliffordConjugation(GA5_4_1 a)	
//(a, -b,-c,-d,-e,-f, -g,-h,-j,-k,-l, -m,-n,-p,-r,-s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, -A)
{
	GA5_4_1 b;

	b.q =  a.q;

	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t =  -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = -a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = -a.xt;		b.yt = -a.yt;		b.zt = -a.zt;

	b.wxy = a.wxy;		b.wxz = a.wxz;		b.wxt = a.wxt;		b.wyz = a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;

	b.wxyzt = -a.wxyzt;

	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 Dual(GA5_4_1 a)	// dual =  a*I_inv 
{

//	GA5_4_1 I_inv,b;	// I_inv created initially zero
//	I_inv.wxyzt = -1;	// pure pseudoscalar
//	b = Product(a,I_inv);	// form dual

	GA5_4_1 b;

//Dual(MV1)  = (A, B,-C,D,-E,-F, -G,H,-J,-K,-L,  M,N,-P,-R,S, -s,r,p,-n,-m,  l,k,j,-h,g, f,e,-d,c,-b, -a)

	b.q = a.wxyzt;

	b.w =  a.xyzt; 	b.x = -a.wyzt; 	b.y =  a.wxzt; 	b.z = -a.wxyt; 	b.t = -a.wxyz;

	b.wx = -a.yzt; 	b.wy =  a.xzt; 	b.wz = -a.xyt; 	b.wt = -a.xyz; 	b.xy = -a.wzt;
	b.xz =  a.wyt; 	b.xt =  a.wyz; 	b.yz = -a.wxt; 	b.yt = -a.wxz; 	b.zt =  a.wxy;

	b.wxy = -a.zt; 	b.wxz =  a.yt; 	b.wxt =  a.yz; 	b.wyz = -a.xt; 	b.wyt = -a.xz;
	b.wzt =  a.xy; 	b.xyz =  a.wt; 	b.xyt =  a.wz; 	b.xzt = -a.wy; 	b.yzt =  a.wx; 

	b.wxyz =  a.t; 	b.wxyt =  a.z; 	b.wxzt = -a.y; 	b.wyzt =  a.x; 	b.xyzt = -a.w;

	b.wxyzt =-a.q;


	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 DorstDual(GA5_4_1 a)
{
	GA5_4_1 b, I_inv;

	I_inv.wxyzt =-1;

	b = LeftContraction(a,I_inv);
	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 DorstUnDual(GA5_4_1 a)
{
	GA5_4_1 b, I;

	I.wxyzt = 1;

	b = LeftContraction(a,I);
	return b;
}


//////////////////////////////////////////////////////

GA5_4_1 Add( GA5_4_1 u,  GA5_4_1 v)
{
	GA5_4_1 w;
	w.q =   u.q   + v.q  ;

	w.w   = u.w   + v.w  ;
	w.x   = u.x   + v.x  ;
	w.y   = u.y   + v.y  ;
	w.z   = u.z   + v.z  ;
	w.t   = u.t   + v.t  ;

	w.wx  = u.wx  + v.wx ;
	w.wy  = u.wy  + v.wy ;
	w.wz  = u.wz  + v.wz ;
	w.wt  = u.wt  + v.wt ;
	w.xy  = u.xy  + v.xy ;
	w.xz  = u.xz  + v.xz ;
	w.yz  = u.yz  + v.yz ;
	w.xt  = u.xt  + v.xt ;
	w.yt  = u.yt  + v.yt ;
	w.zt  = u.zt  + v.zt ;

	w.wxy = u.wxy + v.wxy;
	w.wxz = u.wxz + v.wxz;
	w.wxt = u.wxt + v.wxt;
	w.wyz = u.wyz + v.wyz;
	w.wyt = u.wyt + v.wyt;
	w.wzt = u.wzt + v.wzt;
	w.xyz = u.xyz + v.xyz;
	w.xyt = u.xyt + v.xyt;
	w.xzt = u.xzt + v.xzt;
	w.yzt = u.yzt + v.yzt;

	w.wxyz = u.wxyz + v.wxyz;
	w.wxyt = u.wxyt + v.wxyt;
	w.wxzt = u.wxzt + v.wxzt;
	w.wyzt = u.wyzt + v.wyzt;
	w.xyzt = u.xyzt + v.xyzt;

	w.wxyzt = u.wxyzt + v.wxyzt;

	return w;
}

//////////////////////////////////////////////////////

GA5_4_1 Subtract(GA5_4_1 u,  GA5_4_1 v)
{
	GA5_4_1 w;
	w.q =   u.q   - v.q  ;

	w.w   = u.w   - v.w  ;
	w.x   = u.x   - v.x  ;
	w.y   = u.y   - v.y  ;
	w.z   = u.z   - v.z  ;
	w.t   = u.t   - v.t  ;

	w.wx  = u.wx  - v.wx ;
	w.wy  = u.wy  - v.wy ;
	w.wz  = u.wz  - v.wz ;
	w.wt  = u.wt  - v.wt ;
	w.xy  = u.xy  - v.xy ;
	w.xz  = u.xz  - v.xz ;
	w.yz  = u.yz  - v.yz ;
	w.xt  = u.xt  - v.xt ;
	w.yt  = u.yt  - v.yt ;
	w.zt  = u.zt  - v.zt ;

	w.wxy = u.wxy - v.wxy;
	w.wxz = u.wxz - v.wxz;
	w.wxt = u.wxt - v.wxt;
	w.wyz = u.wyz - v.wyz;
	w.wyt = u.wyt - v.wyt;
	w.wzt = u.wzt - v.wzt;
	w.xyz = u.xyz - v.xyz;
	w.xyt = u.xyt - v.xyt;
	w.xzt = u.xzt - v.xzt;
	w.yzt = u.yzt - v.yzt;

	w.wxyz = u.wxyz - v.wxyz;
	w.wxyt = u.wxyt - v.wxyt;
	w.wxzt = u.wxzt - v.wxzt;
	w.wyzt = u.wyzt - v.wyzt;
	w.xyzt = u.xyzt - v.xyzt;

	w.wxyzt = u.wxyzt - v.wxyzt;

	return w;
}

//////////////////////////////////////////////////////

int Equal(GA5_4_1 u,  GA5_4_1 v)
{
	int result;
	result = 	(u.q ==v.q )&&

			(u.w ==v.w )&&
			(u.x ==v.x )&&
			(u.y ==v.y )&&
			(u.z ==v.z )&&
			(u.t ==v.t )&&

			(u.wx==v.wx)&&
			(u.wy==v.wy)&&
			(u.wz==v.wz)&&
			(u.wt==v.wt)&&
			(u.xy==v.xy)&&
			(u.xz==v.xz)&&
			(u.yz==v.yz)&&
			(u.xt==v.xt)&&
			(u.yt==v.yt)&&
			(u.zt==v.zt)&&

			(u.wxy==v.wxy)&&
			(u.wxz==v.wxz)&&
			(u.wxt==v.wxt)&&
			(u.wyz==v.wyz)&&
			(u.wyt==v.wyt)&&
			(u.wzt==v.wzt)&&
			(u.xyz==v.xyz)&&
			(u.xyt==v.xyt)&&
			(u.xzt==v.xzt)&&
			(u.yzt==v.yzt)&&

			(u.wxyz==v.wxyz)&&
			(u.wxyt==v.wxyt)&&
			(u.wxzt==v.wxzt)&&
			(u.wyzt==v.wyzt)&&
			(u.xyzt==v.xyzt)&&

			(u.wxyzt==v.wxyzt);
	return result;
}

//////////////////////////////////////////////////////

int Not_Equal(GA5_4_1 u,  GA5_4_1 v)
{
	int result;

	result = 	(u.q !=v.q )||

			(u.w !=v.w )||
			(u.x !=v.x )||
			(u.y !=v.y )||
			(u.z !=v.z )||
			(u.t !=v.t )||

			(u.wx!=v.wx)||
			(u.wy!=v.wy)||
			(u.wz!=v.wz)||
			(u.wt!=v.wt)||
			(u.xy!=v.xy)||
			(u.xz!=v.xz)||
			(u.yz!=v.yz)||
			(u.xt!=v.xt)||
			(u.yt!=v.yt)||
			(u.zt!=v.zt)||

			(u.wxy!=v.wxy)||
			(u.wxz!=v.wxz)||
			(u.wxt!=v.wxt)||
			(u.wyz!=v.wyz)||
			(u.wyt!=v.wyt)||
			(u.wzt!=v.wzt)||
			(u.xyz!=v.xyz)||
			(u.xyt!=v.xyt)||
			(u.xzt!=v.xzt)||
			(u.yzt!=v.yzt)||

			(u.wxyz!=v.wxyz)||
			(u.wxyt!=v.wxyt)||
			(u.wxzt!=v.wxzt)||
			(u.wyzt!=v.wyzt)||
			(u.xyzt!=v.xyzt)||

			(u.wxyzt!=v.wxyzt);
	return result;
}

//////////////////////////////////////////////////////

GA5_4_1 Product( GA5_4_1 a,  GA5_4_1 b) {
	GA5_4_1 c;

c.q = a.q*b.q - a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.q*b.w + a.t*b.wt + a.w*b.q - a.wt*b.t + a.wx*b.x + a.wxt*b.xt - a.wxy*b.xy + a.wxyt*b.xyt - a.wxyz*b.xyz - a.wxyzt*b.xyzt - a.wxz*b.xz + a.wxzt*b.xzt + a.wy*b.y + a.wyt*b.yt - a.wyz*b.yz + a.wyzt*b.yzt + a.wz*b.z + a.wzt*b.zt - a.x*b.wx + a.xt*b.wxt - a.xy*b.wxy - a.xyt*b.wxyt + a.xyz*b.wxyz - a.xyzt*b.wxyzt - a.xz*b.wxz - a.xzt*b.wxzt - a.y*b.wy + a.yt*b.wyt - a.yz*b.wyz - a.yzt*b.wyzt - a.z*b.wz + a.zt*b.wzt) ; 

c.x = + (a.q*b.x + a.t*b.xt + a.w*b.wx - a.wt*b.wxt - a.wx*b.w - a.wxt*b.wt + a.wxy*b.wy - a.wxyt*b.wyt + a.wxyz*b.wyz + a.wxyzt*b.wyzt + a.wxz*b.wz - a.wxzt*b.wzt + a.wy*b.wxy + a.wyt*b.wxyt - a.wyz*b.wxyz + a.wyzt*b.wxyzt + a.wz*b.wxz + a.wzt*b.wxzt + a.x*b.q - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt) ; 

c.y = + (a.q*b.y + a.t*b.yt + a.w*b.wy - a.wt*b.wyt - a.wx*b.wxy - a.wxt*b.wxyt - a.wxy*b.wx + a.wxyt*b.wxt - a.wxyz*b.wxz - a.wxyzt*b.wxzt + a.wxz*b.wxyz - a.wxzt*b.wxyzt - a.wy*b.w - a.wyt*b.wt + a.wyz*b.wz - a.wyzt*b.wzt + a.wz*b.wyz + a.wzt*b.wyzt + a.x*b.xy - a.xt*b.xyt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.xz*b.xyz + a.xzt*b.xyzt + a.y*b.q - a.yt*b.t + a.yz*b.z + a.yzt*b.zt - a.z*b.yz + a.zt*b.yzt) ; 

c.z = + (a.q*b.z + a.t*b.zt + a.w*b.wz - a.wt*b.wzt - a.wx*b.wxz - a.wxt*b.wxzt - a.wxy*b.wxyz + a.wxyt*b.wxyzt + a.wxyz*b.wxy + a.wxyzt*b.wxyt - a.wxz*b.wx + a.wxzt*b.wxt - a.wy*b.wyz - a.wyt*b.wyzt - a.wyz*b.wy + a.wyzt*b.wyt - a.wz*b.w - a.wzt*b.wt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt + a.y*b.yz - a.yt*b.yzt - a.yz*b.y - a.yzt*b.yt + a.z*b.q - a.zt*b.t) ; 

c.t = + (a.q*b.t + a.t*b.q + a.w*b.wt - a.wt*b.w - a.wx*b.wxt - a.wxt*b.wx - a.wxy*b.wxyt + a.wxyt*b.wxy + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz - a.wxz*b.wxzt + a.wxzt*b.wxz - a.wy*b.wyt - a.wyt*b.wy - a.wyz*b.wyzt + a.wyzt*b.wyz - a.wz*b.wzt - a.wzt*b.wz + a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z) ; 



c.wx = + (a.q*b.wx - a.t*b.wxt + a.w*b.x + a.wt*b.xt + a.wx*b.q - a.wxt*b.t + a.wxy*b.y + a.wxyt*b.yt - a.wxyz*b.yz + a.wxyzt*b.yzt + a.wxz*b.z + a.wxzt*b.zt - a.wy*b.xy + a.wyt*b.xyt - a.wyz*b.xyz - a.wyzt*b.xyzt - a.wz*b.xz + a.wzt*b.xzt - a.x*b.w - a.xt*b.wt + a.xy*b.wy - a.xyt*b.wyt + a.xyz*b.wyz + a.xyzt*b.wyzt + a.xz*b.wz - a.xzt*b.wzt + a.y*b.wxy + a.yt*b.wxyt - a.yz*b.wxyz + a.yzt*b.wxyzt + a.z*b.wxz + a.zt*b.wxzt) ; 

c.wy = + (a.q*b.wy - a.t*b.wyt + a.w*b.y + a.wt*b.yt + a.wx*b.xy - a.wxt*b.xyt - a.wxy*b.x - a.wxyt*b.xt + a.wxyz*b.xz - a.wxyzt*b.xzt + a.wxz*b.xyz + a.wxzt*b.xyzt + a.wy*b.q - a.wyt*b.t + a.wyz*b.z + a.wyzt*b.zt - a.wz*b.yz + a.wzt*b.yzt - a.x*b.wxy - a.xt*b.wxyt - a.xy*b.wx + a.xyt*b.wxt - a.xyz*b.wxz - a.xyzt*b.wxzt + a.xz*b.wxyz - a.xzt*b.wxyzt - a.y*b.w - a.yt*b.wt + a.yz*b.wz - a.yzt*b.wzt + a.z*b.wyz + a.zt*b.wyzt) ; 

c.wz = + (a.q*b.wz - a.t*b.wzt + a.w*b.z + a.wt*b.zt + a.wx*b.xz - a.wxt*b.xzt - a.wxy*b.xyz - a.wxyt*b.xyzt - a.wxyz*b.xy + a.wxyzt*b.xyt - a.wxz*b.x - a.wxzt*b.xt + a.wy*b.yz - a.wyt*b.yzt - a.wyz*b.y - a.wyzt*b.yt + a.wz*b.q - a.wzt*b.t - a.x*b.wxz - a.xt*b.wxzt - a.xy*b.wxyz + a.xyt*b.wxyzt + a.xyz*b.wxy + a.xyzt*b.wxyt - a.xz*b.wx + a.xzt*b.wxt - a.y*b.wyz - a.yt*b.wyzt - a.yz*b.wy + a.yzt*b.wyt - a.z*b.w - a.zt*b.wt) ; 

c.wt = + (a.q*b.wt - a.t*b.w + a.w*b.t + a.wt*b.q + a.wx*b.xt - a.wxt*b.x - a.wxy*b.xyt - a.wxyt*b.xy - a.wxyz*b.xyzt + a.wxyzt*b.xyz - a.wxz*b.xzt - a.wxzt*b.xz + a.wy*b.yt - a.wyt*b.y - a.wyz*b.yzt - a.wyzt*b.yz + a.wz*b.zt - a.wzt*b.z - a.x*b.wxt - a.xt*b.wx - a.xy*b.wxyt + a.xyt*b.wxy + a.xyz*b.wxyzt + a.xyzt*b.wxyz - a.xz*b.wxzt + a.xzt*b.wxz - a.y*b.wyt - a.yt*b.wy - a.yz*b.wyzt + a.yzt*b.wyz - a.z*b.wzt - a.zt*b.wz) ; 

c.xy = + (a.q*b.xy - a.t*b.xyt + a.w*b.wxy + a.wt*b.wxyt - a.wx*b.wy + a.wxt*b.wyt + a.wxy*b.w + a.wxyt*b.wt - a.wxyz*b.wz + a.wxyzt*b.wzt - a.wxz*b.wyz - a.wxzt*b.wyzt + a.wy*b.wx - a.wyt*b.wxt + a.wyz*b.wxz + a.wyzt*b.wxzt - a.wz*b.wxyz + a.wzt*b.wxyzt + a.x*b.y + a.xt*b.yt + a.xy*b.q - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt - a.xz*b.yz + a.xzt*b.yzt - a.y*b.x - a.yt*b.xt + a.yz*b.xz - a.yzt*b.xzt + a.z*b.xyz + a.zt*b.xyzt) ; 

c.xz = + (a.q*b.xz - a.t*b.xzt + a.w*b.wxz + a.wt*b.wxzt - a.wx*b.wz + a.wxt*b.wzt + a.wxy*b.wyz + a.wxyt*b.wyzt + a.wxyz*b.wy - a.wxyzt*b.wyt + a.wxz*b.w + a.wxzt*b.wt + a.wy*b.wxyz - a.wyt*b.wxyzt - a.wyz*b.wxy - a.wyzt*b.wxyt + a.wz*b.wx - a.wzt*b.wxt + a.x*b.z + a.xt*b.zt + a.xy*b.yz - a.xyt*b.yzt - a.xyz*b.y - a.xyzt*b.yt + a.xz*b.q - a.xzt*b.t - a.y*b.xyz - a.yt*b.xyzt - a.yz*b.xy + a.yzt*b.xyt - a.z*b.x - a.zt*b.xt) ; 

c.xt = + (a.q*b.xt - a.t*b.x + a.w*b.wxt + a.wt*b.wx - a.wx*b.wt + a.wxt*b.w + a.wxy*b.wyt + a.wxyt*b.wy + a.wxyz*b.wyzt - a.wxyzt*b.wyz + a.wxz*b.wzt + a.wxzt*b.wz + a.wy*b.wxyt - a.wyt*b.wxy - a.wyz*b.wxyzt - a.wyzt*b.wxyz + a.wz*b.wxzt - a.wzt*b.wxz + a.x*b.t + a.xt*b.q + a.xy*b.yt - a.xyt*b.y - a.xyz*b.yzt - a.xyzt*b.yz + a.xz*b.zt - a.xzt*b.z - a.y*b.xyt - a.yt*b.xy - a.yz*b.xyzt + a.yzt*b.xyz - a.z*b.xzt - a.zt*b.xz) ; 

c.yz = + (a.q*b.yz - a.t*b.yzt + a.w*b.wyz + a.wt*b.wyzt - a.wx*b.wxyz + a.wxt*b.wxyzt - a.wxy*b.wxz - a.wxyt*b.wxzt - a.wxyz*b.wx + a.wxyzt*b.wxt + a.wxz*b.wxy + a.wxzt*b.wxyt - a.wy*b.wz + a.wyt*b.wzt + a.wyz*b.w + a.wyzt*b.wt + a.wz*b.wy - a.wzt*b.wyt + a.x*b.xyz + a.xt*b.xyzt - a.xy*b.xz + a.xyt*b.xzt + a.xyz*b.x + a.xyzt*b.xt + a.xz*b.xy - a.xzt*b.xyt + a.y*b.z + a.yt*b.zt + a.yz*b.q - a.yzt*b.t - a.z*b.y - a.zt*b.yt) ; 

c.yt = + (a.q*b.yt - a.t*b.y + a.w*b.wyt + a.wt*b.wy - a.wx*b.wxyt + a.wxt*b.wxy - a.wxy*b.wxt - a.wxyt*b.wx - a.wxyz*b.wxzt + a.wxyzt*b.wxz + a.wxz*b.wxyzt + a.wxzt*b.wxyz - a.wy*b.wt + a.wyt*b.w + a.wyz*b.wzt + a.wyzt*b.wz + a.wz*b.wyzt - a.wzt*b.wyz + a.x*b.xyt + a.xt*b.xy - a.xy*b.xt + a.xyt*b.x + a.xyz*b.xzt + a.xyzt*b.xz + a.xz*b.xyzt - a.xzt*b.xyz + a.y*b.t + a.yt*b.q + a.yz*b.zt - a.yzt*b.z - a.z*b.yzt - a.zt*b.yz) ; 

c.zt = + (a.q*b.zt - a.t*b.z + a.w*b.wzt + a.wt*b.wz - a.wx*b.wxzt + a.wxt*b.wxz - a.wxy*b.wxyzt - a.wxyt*b.wxyz + a.wxyz*b.wxyt - a.wxyzt*b.wxy - a.wxz*b.wxt - a.wxzt*b.wx - a.wy*b.wyzt + a.wyt*b.wyz - a.wyz*b.wyt - a.wyzt*b.wy - a.wz*b.wt + a.wzt*b.w + a.x*b.xzt + a.xt*b.xz - a.xy*b.xyzt + a.xyt*b.xyz - a.xyz*b.xyt - a.xyzt*b.xy - a.xz*b.xt + a.xzt*b.x + a.y*b.yzt + a.yt*b.yz - a.yz*b.yt + a.yzt*b.y + a.z*b.t + a.zt*b.q) ; 



c.wxy = + (a.q*b.wxy + a.t*b.wxyt + a.w*b.xy - a.wt*b.xyt + a.wx*b.y + a.wxt*b.yt + a.wxy*b.q - a.wxyt*b.t + a.wxyz*b.z + a.wxyzt*b.zt - a.wxz*b.yz + a.wxzt*b.yzt - a.wy*b.x - a.wyt*b.xt + a.wyz*b.xz - a.wyzt*b.xzt + a.wz*b.xyz + a.wzt*b.xyzt - a.x*b.wy + a.xt*b.wyt + a.xy*b.w + a.xyt*b.wt - a.xyz*b.wz + a.xyzt*b.wzt - a.xz*b.wyz - a.xzt*b.wyzt + a.y*b.wx - a.yt*b.wxt + a.yz*b.wxz + a.yzt*b.wxzt - a.z*b.wxyz + a.zt*b.wxyzt) ; 

c.wxz = + (a.q*b.wxz + a.t*b.wxzt + a.w*b.xz - a.wt*b.xzt + a.wx*b.z + a.wxt*b.zt + a.wxy*b.yz - a.wxyt*b.yzt - a.wxyz*b.y - a.wxyzt*b.yt + a.wxz*b.q - a.wxzt*b.t - a.wy*b.xyz - a.wyt*b.xyzt - a.wyz*b.xy + a.wyzt*b.xyt - a.wz*b.x - a.wzt*b.xt - a.x*b.wz + a.xt*b.wzt + a.xy*b.wyz + a.xyt*b.wyzt + a.xyz*b.wy - a.xyzt*b.wyt + a.xz*b.w + a.xzt*b.wt + a.y*b.wxyz - a.yt*b.wxyzt - a.yz*b.wxy - a.yzt*b.wxyt + a.z*b.wx - a.zt*b.wxt) ; 

c.wxt = + (a.q*b.wxt + a.t*b.wx + a.w*b.xt - a.wt*b.x + a.wx*b.t + a.wxt*b.q + a.wxy*b.yt - a.wxyt*b.y - a.wxyz*b.yzt - a.wxyzt*b.yz + a.wxz*b.zt - a.wxzt*b.z - a.wy*b.xyt - a.wyt*b.xy - a.wyz*b.xyzt + a.wyzt*b.xyz - a.wz*b.xzt - a.wzt*b.xz - a.x*b.wt + a.xt*b.w + a.xy*b.wyt + a.xyt*b.wy + a.xyz*b.wyzt - a.xyzt*b.wyz + a.xz*b.wzt + a.xzt*b.wz + a.y*b.wxyt - a.yt*b.wxy - a.yz*b.wxyzt - a.yzt*b.wxyz + a.z*b.wxzt - a.zt*b.wxz) ; 

c.wyz = + (a.q*b.wyz + a.t*b.wyzt + a.w*b.yz - a.wt*b.yzt + a.wx*b.xyz + a.wxt*b.xyzt - a.wxy*b.xz + a.wxyt*b.xzt + a.wxyz*b.x + a.wxyzt*b.xt + a.wxz*b.xy - a.wxzt*b.xyt + a.wy*b.z + a.wyt*b.zt + a.wyz*b.q - a.wyzt*b.t - a.wz*b.y - a.wzt*b.yt - a.x*b.wxyz + a.xt*b.wxyzt - a.xy*b.wxz - a.xyt*b.wxzt - a.xyz*b.wx + a.xyzt*b.wxt + a.xz*b.wxy + a.xzt*b.wxyt - a.y*b.wz + a.yt*b.wzt + a.yz*b.w + a.yzt*b.wt + a.z*b.wy - a.zt*b.wyt) ; 

c.wyt = + (a.q*b.wyt + a.t*b.wy + a.w*b.yt - a.wt*b.y + a.wx*b.xyt + a.wxt*b.xy - a.wxy*b.xt + a.wxyt*b.x + a.wxyz*b.xzt + a.wxyzt*b.xz + a.wxz*b.xyzt - a.wxzt*b.xyz + a.wy*b.t + a.wyt*b.q + a.wyz*b.zt - a.wyzt*b.z - a.wz*b.yzt - a.wzt*b.yz - a.x*b.wxyt + a.xt*b.wxy - a.xy*b.wxt - a.xyt*b.wx - a.xyz*b.wxzt + a.xyzt*b.wxz + a.xz*b.wxyzt + a.xzt*b.wxyz - a.y*b.wt + a.yt*b.w + a.yz*b.wzt + a.yzt*b.wz + a.z*b.wyzt - a.zt*b.wyz) ; 

c.wzt = + (a.q*b.wzt + a.t*b.wz + a.w*b.zt - a.wt*b.z + a.wx*b.xzt + a.wxt*b.xz - a.wxy*b.xyzt + a.wxyt*b.xyz - a.wxyz*b.xyt - a.wxyzt*b.xy - a.wxz*b.xt + a.wxzt*b.x + a.wy*b.yzt + a.wyt*b.yz - a.wyz*b.yt + a.wyzt*b.y + a.wz*b.t + a.wzt*b.q - a.x*b.wxzt + a.xt*b.wxz - a.xy*b.wxyzt - a.xyt*b.wxyz + a.xyz*b.wxyt - a.xyzt*b.wxy - a.xz*b.wxt - a.xzt*b.wx - a.y*b.wyzt + a.yt*b.wyz - a.yz*b.wyt - a.yzt*b.wy - a.z*b.wt + a.zt*b.w) ; 

c.xyz = + (a.q*b.xyz + a.t*b.xyzt + a.w*b.wxyz - a.wt*b.wxyzt - a.wx*b.wyz - a.wxt*b.wyzt + a.wxy*b.wz - a.wxyt*b.wzt - a.wxyz*b.w - a.wxyzt*b.wt - a.wxz*b.wy + a.wxzt*b.wyt + a.wy*b.wxz + a.wyt*b.wxzt + a.wyz*b.wx - a.wyzt*b.wxt - a.wz*b.wxy - a.wzt*b.wxyt + a.x*b.yz - a.xt*b.yzt + a.xy*b.z + a.xyt*b.zt + a.xyz*b.q - a.xyzt*b.t - a.xz*b.y - a.xzt*b.yt - a.y*b.xz + a.yt*b.xzt + a.yz*b.x + a.yzt*b.xt + a.z*b.xy - a.zt*b.xyt) ; 

c.xyt = + (a.q*b.xyt + a.t*b.xy + a.w*b.wxyt - a.wt*b.wxy - a.wx*b.wyt - a.wxt*b.wy + a.wxy*b.wt - a.wxyt*b.w - a.wxyz*b.wzt - a.wxyzt*b.wz - a.wxz*b.wyzt + a.wxzt*b.wyz + a.wy*b.wxt + a.wyt*b.wx + a.wyz*b.wxzt - a.wyzt*b.wxz - a.wz*b.wxyzt - a.wzt*b.wxyz + a.x*b.yt - a.xt*b.y + a.xy*b.t + a.xyt*b.q + a.xyz*b.zt - a.xyzt*b.z - a.xz*b.yzt - a.xzt*b.yz - a.y*b.xt + a.yt*b.x + a.yz*b.xzt + a.yzt*b.xz + a.z*b.xyzt - a.zt*b.xyz) ; 

c.xzt = + (a.q*b.xzt + a.t*b.xz + a.w*b.wxzt - a.wt*b.wxz - a.wx*b.wzt - a.wxt*b.wz + a.wxy*b.wyzt - a.wxyt*b.wyz + a.wxyz*b.wyt + a.wxyzt*b.wy + a.wxz*b.wt - a.wxzt*b.w + a.wy*b.wxyzt + a.wyt*b.wxyz - a.wyz*b.wxyt + a.wyzt*b.wxy + a.wz*b.wxt + a.wzt*b.wx + a.x*b.zt - a.xt*b.z + a.xy*b.yzt + a.xyt*b.yz - a.xyz*b.yt + a.xyzt*b.y + a.xz*b.t + a.xzt*b.q - a.y*b.xyzt + a.yt*b.xyz - a.yz*b.xyt - a.yzt*b.xy - a.z*b.xt + a.zt*b.x) ; 

c.yzt = + (a.q*b.yzt + a.t*b.yz + a.w*b.wyzt - a.wt*b.wyz - a.wx*b.wxyzt - a.wxt*b.wxyz - a.wxy*b.wxzt + a.wxyt*b.wxz - a.wxyz*b.wxt - a.wxyzt*b.wx + a.wxz*b.wxyt - a.wxzt*b.wxy - a.wy*b.wzt - a.wyt*b.wz + a.wyz*b.wt - a.wyzt*b.w + a.wz*b.wyt + a.wzt*b.wy + a.x*b.xyzt - a.xt*b.xyz - a.xy*b.xzt - a.xyt*b.xz + a.xyz*b.xt - a.xyzt*b.x + a.xz*b.xyt + a.xzt*b.xy + a.y*b.zt - a.yt*b.z + a.yz*b.t + a.yzt*b.q - a.z*b.yt + a.zt*b.y) ; 



c.wxyz = + (a.q*b.wxyz - a.t*b.wxyzt + a.w*b.xyz + a.wt*b.xyzt + a.wx*b.yz - a.wxt*b.yzt + a.wxy*b.z + a.wxyt*b.zt + a.wxyz*b.q - a.wxyzt*b.t - a.wxz*b.y - a.wxzt*b.yt - a.wy*b.xz + a.wyt*b.xzt + a.wyz*b.x + a.wyzt*b.xt + a.wz*b.xy - a.wzt*b.xyt - a.x*b.wyz - a.xt*b.wyzt + a.xy*b.wz - a.xyt*b.wzt - a.xyz*b.w - a.xyzt*b.wt - a.xz*b.wy + a.xzt*b.wyt + a.y*b.wxz + a.yt*b.wxzt + a.yz*b.wx - a.yzt*b.wxt - a.z*b.wxy - a.zt*b.wxyt) ; 

c.wxyt = + (a.q*b.wxyt - a.t*b.wxy + a.w*b.xyt + a.wt*b.xy + a.wx*b.yt - a.wxt*b.y + a.wxy*b.t + a.wxyt*b.q + a.wxyz*b.zt - a.wxyzt*b.z - a.wxz*b.yzt - a.wxzt*b.yz - a.wy*b.xt + a.wyt*b.x + a.wyz*b.xzt + a.wyzt*b.xz + a.wz*b.xyzt - a.wzt*b.xyz - a.x*b.wyt - a.xt*b.wy + a.xy*b.wt - a.xyt*b.w - a.xyz*b.wzt - a.xyzt*b.wz - a.xz*b.wyzt + a.xzt*b.wyz + a.y*b.wxt + a.yt*b.wx + a.yz*b.wxzt - a.yzt*b.wxz - a.z*b.wxyzt - a.zt*b.wxyz) ; 

c.wxzt = + (a.q*b.wxzt - a.t*b.wxz + a.w*b.xzt + a.wt*b.xz + a.wx*b.zt - a.wxt*b.z + a.wxy*b.yzt + a.wxyt*b.yz - a.wxyz*b.yt + a.wxyzt*b.y + a.wxz*b.t + a.wxzt*b.q - a.wy*b.xyzt + a.wyt*b.xyz - a.wyz*b.xyt - a.wyzt*b.xy - a.wz*b.xt + a.wzt*b.x - a.x*b.wzt - a.xt*b.wz + a.xy*b.wyzt - a.xyt*b.wyz + a.xyz*b.wyt + a.xyzt*b.wy + a.xz*b.wt - a.xzt*b.w + a.y*b.wxyzt + a.yt*b.wxyz - a.yz*b.wxyt + a.yzt*b.wxy + a.z*b.wxt + a.zt*b.wx) ; 

c.wyzt = + (a.q*b.wyzt - a.t*b.wyz + a.w*b.yzt + a.wt*b.yz + a.wx*b.xyzt - a.wxt*b.xyz - a.wxy*b.xzt - a.wxyt*b.xz + a.wxyz*b.xt - a.wxyzt*b.x + a.wxz*b.xyt + a.wxzt*b.xy + a.wy*b.zt - a.wyt*b.z + a.wyz*b.t + a.wyzt*b.q - a.wz*b.yt + a.wzt*b.y - a.x*b.wxyzt - a.xt*b.wxyz - a.xy*b.wxzt + a.xyt*b.wxz - a.xyz*b.wxt - a.xyzt*b.wx + a.xz*b.wxyt - a.xzt*b.wxy - a.y*b.wzt - a.yt*b.wz + a.yz*b.wt - a.yzt*b.w + a.z*b.wyt + a.zt*b.wy) ; 

c.xyzt = + (a.q*b.xyzt - a.t*b.xyz + a.w*b.wxyzt + a.wt*b.wxyz - a.wx*b.wyzt + a.wxt*b.wyz + a.wxy*b.wzt + a.wxyt*b.wz - a.wxyz*b.wt + a.wxyzt*b.w - a.wxz*b.wyt - a.wxzt*b.wy + a.wy*b.wxzt - a.wyt*b.wxz + a.wyz*b.wxt + a.wyzt*b.wx - a.wz*b.wxyt + a.wzt*b.wxy + a.x*b.yzt + a.xt*b.yz + a.xy*b.zt - a.xyt*b.z + a.xyz*b.t + a.xyzt*b.q - a.xz*b.yt + a.xzt*b.y - a.y*b.xzt - a.yt*b.xz + a.yz*b.xt - a.yzt*b.x + a.z*b.xyt + a.zt*b.xy) ; 


c.wxyzt = + (a.q*b.wxyzt + a.t*b.wxyz + a.w*b.xyzt - a.wt*b.xyz + a.wx*b.yzt + a.wxt*b.yz + a.wxy*b.zt - a.wxyt*b.z + a.wxyz*b.t + a.wxyzt*b.q - a.wxz*b.yt + a.wxzt*b.y - a.wy*b.xzt - a.wyt*b.xz + a.wyz*b.xt - a.wyzt*b.x + a.wz*b.xyt + a.wzt*b.xy - a.x*b.wyzt + a.xt*b.wyz + a.xy*b.wzt + a.xyt*b.wz - a.xyz*b.wt + a.xyzt*b.w - a.xz*b.wyt - a.xzt*b.wy + a.y*b.wxzt - a.yt*b.wxz + a.yz*b.wxt + a.yzt*b.wx - a.z*b.wxyt + a.zt*b.wxy) ; 

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Divide_By_Constant(GA5_4_1 u,  long i)
{
	GA5_4_1 a;
	a.q = u.q/i;

	a.w = u.w/i;
	a.x = u.x/i;
	a.y = u.y/i;
	a.z = u.z/i;
	a.t = u.t/i;

	a.wx = u.wx/i;
	a.wy = u.wy/i;
	a.wz = u.wz/i;
	a.wt = u.wt/i;
	a.xy = u.xy/i;
	a.xz = u.xz/i;
	a.yz = u.yz/i;
	a.xt = u.xt/i;
	a.yt = u.yt/i;
	a.zt = u.zt/i;

	a.wxy = u.wxy/i;
	a.wxz = u.wxz/i;
	a.wxt = u.wxt/i;
	a.wyz = u.wyz/i;
	a.wyt = u.wyt/i;
	a.wzt = u.wzt/i;
	a.xyz = u.xyz/i;
	a.xyt = u.xyt/i;
	a.xzt = u.xzt/i;
	a.yzt = u.yzt/i;

	a.wxyz = u.wxyz/i;
	a.wxyt = u.wxyt/i;
	a.wxzt = u.wxzt/i;
	a.wyzt = u.wyzt/i;
	a.xyzt = u.xyzt/i;

	a.wxyzt = u.wxyzt/i;

	return a;
}

//////////////////////////////////////////////////////

GA5_4_1 Wedge(GA5_4_1 a,  GA5_4_1 b) {

	GA5_4_1 c;

c.q    =  + a.q   *b.q    ; 
c.w    =  + a.q   *b.w    + a.w   *b.q    ; 
c.x    =  + a.q   *b.x    + a.x   *b.q    ; 
c.y    =  + a.q   *b.y    + a.y   *b.q    ; 
c.z    =  + a.q   *b.z    + a.z   *b.q    ; 
c.t    =  + a.q   *b.t    + a.t   *b.q    ; 
c.wx   =  + a.q   *b.wx   + a.w   *b.x    - a.x   *b.w    + a.wx  *b.q    ; 
c.wy   =  + a.q   *b.wy   + a.w   *b.y    - a.y   *b.w    + a.wy  *b.q    ; 
c.wz   =  + a.q   *b.wz   + a.w   *b.z    - a.z   *b.w    + a.wz  *b.q    ; 
c.wt   =  + a.q   *b.wt   + a.w   *b.t    - a.t   *b.w    + a.wt  *b.q    ; 
c.xy   =  + a.q   *b.xy   + a.x   *b.y    - a.y   *b.x    + a.xy  *b.q    ; 
c.xz   =  + a.q   *b.xz   + a.x   *b.z    - a.z   *b.x    + a.xz  *b.q    ; 
c.xt   =  + a.q   *b.xt   + a.x   *b.t    - a.t   *b.x    + a.xt  *b.q    ; 
c.yz   =  + a.q   *b.yz   + a.y   *b.z    - a.z   *b.y    + a.yz  *b.q    ; 
c.yt   =  + a.q   *b.yt   + a.y   *b.t    - a.t   *b.y    + a.yt  *b.q    ; 
c.zt   =  + a.q   *b.zt   + a.z   *b.t    - a.t   *b.z    + a.zt  *b.q    ; 
c.wxy  =  + a.q   *b.wxy  + a.w   *b.xy   - a.x   *b.wy   + a.y   *b.wx   + a.wx  *b.y    - a.wy  *b.x    + a.xy  *b.w    + a.wxy *b.q    ; 
c.wxz  =  + a.q   *b.wxz  + a.w   *b.xz   - a.x   *b.wz   + a.z   *b.wx   + a.wx  *b.z    - a.wz  *b.x    + a.xz  *b.w    + a.wxz *b.q    ; 
c.wxt  =  + a.q   *b.wxt  + a.w   *b.xt   - a.x   *b.wt   + a.t   *b.wx   + a.wx  *b.t    - a.wt  *b.x    + a.xt  *b.w    + a.wxt *b.q    ; 
c.wyz  =  + a.q   *b.wyz  + a.w   *b.yz   - a.y   *b.wz   + a.z   *b.wy   + a.wy  *b.z    - a.wz  *b.y    + a.yz  *b.w    + a.wyz *b.q    ; 
c.wyt  =  + a.q   *b.wyt  + a.w   *b.yt   - a.y   *b.wt   + a.t   *b.wy   + a.wy  *b.t    - a.wt  *b.y    + a.yt  *b.w    + a.wyt *b.q    ; 
c.wzt  =  + a.q   *b.wzt  + a.w   *b.zt   - a.z   *b.wt   + a.t   *b.wz   + a.wz  *b.t    - a.wt  *b.z    + a.zt  *b.w    + a.wzt *b.q    ; 
c.xyz  =  + a.q   *b.xyz  + a.x   *b.yz   - a.y   *b.xz   + a.z   *b.xy   + a.xy  *b.z    - a.xz  *b.y    + a.yz  *b.x    + a.xyz *b.q    ; 
c.xyt  =  + a.q   *b.xyt  + a.x   *b.yt   - a.y   *b.xt   + a.t   *b.xy   + a.xy  *b.t    - a.xt  *b.y    + a.yt  *b.x    + a.xyt *b.q    ; 
c.xzt  =  + a.q   *b.xzt  + a.x   *b.zt   - a.z   *b.xt   + a.t   *b.xz   + a.xz  *b.t    - a.xt  *b.z    + a.zt  *b.x    + a.xzt *b.q    ; 
c.yzt  =  + a.q   *b.yzt  + a.y   *b.zt   - a.z   *b.yt   + a.t   *b.yz   + a.yz  *b.t    - a.yt  *b.z    + a.zt  *b.y    + a.yzt *b.q    ; 
c.wxyz =  + a.q   *b.wxyz + a.w   *b.xyz  - a.x   *b.wyz  + a.y   *b.wxz  - a.z   *b.wxy  + a.wx  *b.yz   - a.wy  *b.xz   + a.wz  *b.xy   + a.xy  *b.wz   - a.xz  *b.wy   + a.yz  *b.wx   + a.wxy *b.z    - a.wxz *b.y    + a.wyz *b.x    - a.xyz *b.w    + a.wxyz*b.q    ; 
c.wxyt =  + a.q   *b.wxyt + a.w   *b.xyt  - a.x   *b.wyt  + a.y   *b.wxt  - a.t   *b.wxy  + a.wx  *b.yt   - a.wy  *b.xt   + a.wt  *b.xy   + a.xy  *b.wt   - a.xt  *b.wy   + a.yt  *b.wx   + a.wxy *b.t    - a.wxt *b.y    + a.wyt *b.x    - a.xyt *b.w    + a.wxyt*b.q    ; 
c.wxzt =  + a.q   *b.wxzt + a.w   *b.xzt  - a.x   *b.wzt  + a.z   *b.wxt  - a.t   *b.wxz  + a.wx  *b.zt   - a.wz  *b.xt   + a.wt  *b.xz   + a.xz  *b.wt   - a.xt  *b.wz   + a.zt  *b.wx   + a.wxz *b.t    - a.wxt *b.z    + a.wzt *b.x    - a.xzt *b.w    + a.wxzt*b.q    ; 
c.wyzt =  + a.q   *b.wyzt + a.w   *b.yzt  - a.y   *b.wzt  + a.z   *b.wyt  - a.t   *b.wyz  + a.wy  *b.zt   - a.wz  *b.yt   + a.wt  *b.yz   + a.yz  *b.wt   - a.yt  *b.wz   + a.zt  *b.wy   + a.wyz *b.t    - a.wyt *b.z    + a.wzt *b.y    - a.yzt *b.w    + a.wyzt*b.q    ; 
c.xyzt =  + a.q   *b.xyzt + a.x   *b.yzt  - a.y   *b.xzt  + a.z   *b.xyt  - a.t   *b.xyz  + a.xy  *b.zt   - a.xz  *b.yt   + a.xt  *b.yz   + a.yz  *b.xt   - a.yt  *b.xz   + a.zt  *b.xy   + a.xyz *b.t    - a.xyt *b.z    + a.xzt *b.y    - a.yzt *b.x    + a.xyzt*b.q    ; 
c.wxyzt =  + a.q   *b.wxyzt + a.w   *b.xyzt - a.x   *b.wyzt + a.y   *b.wxzt - a.z   *b.wxyt + a.t   *b.wxyz + a.wx  *b.yzt  - a.wy  *b.xzt  + a.wz  *b.xyt  - a.wt  *b.xyz  + a.xy  *b.wzt  - a.xz  *b.wyt  + a.xt  *b.wyz  + a.yz  *b.wxt  - a.yt  *b.wxz  + a.zt  *b.wxy  + a.wxy *b.zt   - a.wxz *b.yt   + a.wxt *b.yz   + a.wyz *b.xt   - a.wyt *b.xz   + a.wzt *b.xy   - a.xyz *b.wt   + a.xyt *b.wz   - a.xzt *b.wy   + a.yzt *b.wx   + a.wxyz*b.t    - a.wxyt*b.z    + a.wxzt*b.y    - a.wyzt*b.x    + a.xyzt*b.w    + a.wxyzt*b.q    ; 

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 AntiWedge(GA5_4_1 a,  GA5_4_1 b) {

	GA5_4_1 c;

//	c = OverBar(Wedge(UnderBar(a),UnderBar(b)));

	c.q =  + a.q*b.wxyzt + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz
		 + a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt
		 + a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt
		 + a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz
		 - a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x
		 + a.xyzt*b.w + a.wxyzt*b.q ; 
	c.w =  + a.w*b.wxyzt + a.wx*b.wyzt - a.wy*b.wxzt + a.wz*b.wxyt - a.wt*b.wxyz + a.wxy*b.wzt
		 - a.wxz*b.wyt + a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz + a.wzt*b.wxy + a.wxyz*b.wt
		 - a.wxyt*b.wz + a.wxzt*b.wy - a.wyzt*b.wx + a.wxyzt*b.w ; 
	c.x =  + a.x*b.wxyzt + a.wx*b.xyzt - a.xy*b.wxzt + a.xz*b.wxyt - a.xt*b.wxyz + a.wxy*b.xzt
		 - a.wxz*b.xyt + a.wxt*b.xyz + a.xyz*b.wxt - a.xyt*b.wxz + a.xzt*b.wxy + a.wxyz*b.xt
		 - a.wxyt*b.xz + a.wxzt*b.xy - a.xyzt*b.wx + a.wxyzt*b.x ; 
	c.y =  + a.y*b.wxyzt + a.wy*b.xyzt - a.xy*b.wyzt + a.yz*b.wxyt - a.yt*b.wxyz + a.wxy*b.yzt
		 - a.wyz*b.xyt + a.wyt*b.xyz + a.xyz*b.wyt - a.xyt*b.wyz + a.yzt*b.wxy + a.wxyz*b.yt
		 - a.wxyt*b.yz + a.wyzt*b.xy - a.xyzt*b.wy + a.wxyzt*b.y ; 
	c.z =  + a.z*b.wxyzt + a.wz*b.xyzt - a.xz*b.wyzt + a.yz*b.wxzt - a.zt*b.wxyz + a.wxz*b.yzt
		 - a.wyz*b.xzt + a.wzt*b.xyz + a.xyz*b.wzt - a.xzt*b.wyz + a.yzt*b.wxz + a.wxyz*b.zt
		 - a.wxzt*b.yz + a.wyzt*b.xz - a.xyzt*b.wz + a.wxyzt*b.z ; 
	c.t =  + a.t*b.wxyzt + a.wt*b.xyzt - a.xt*b.wyzt + a.yt*b.wxzt - a.zt*b.wxyt + a.wxt*b.yzt
		 - a.wyt*b.xzt + a.wzt*b.xyt + a.xyt*b.wzt - a.xzt*b.wyt + a.yzt*b.wxt + a.wxyt*b.zt
		 - a.wxzt*b.yt + a.wyzt*b.xt - a.xyzt*b.wt + a.wxyzt*b.t ; 
	
	c.wx =  + a.wx*b.wxyzt + a.wxy*b.wxzt - a.wxz*b.wxyt + a.wxt*b.wxyz + a.wxyz*b.wxt - a.wxyt*b.wxz + a.wxzt*b.wxy + a.wxyzt*b.wx ; 
	c.wy =  + a.wy*b.wxyzt + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz + a.wxyz*b.wyt - a.wxyt*b.wyz + a.wyzt*b.wxy + a.wxyzt*b.wy ; 
	c.wz =  + a.wz*b.wxyzt + a.wxz*b.wyzt - a.wyz*b.wxzt + a.wzt*b.wxyz + a.wxyz*b.wzt - a.wxzt*b.wyz + a.wyzt*b.wxz + a.wxyzt*b.wz ; 
	c.wt =  + a.wt*b.wxyzt + a.wxt*b.wyzt - a.wyt*b.wxzt + a.wzt*b.wxyt + a.wxyt*b.wzt - a.wxzt*b.wyt + a.wyzt*b.wxt + a.wxyzt*b.wt ; 
	c.xy =  + a.xy*b.wxyzt + a.wxy*b.xyzt - a.xyz*b.wxyt + a.xyt*b.wxyz + a.wxyz*b.xyt - a.wxyt*b.xyz + a.xyzt*b.wxy + a.wxyzt*b.xy ; 
	c.xz =  + a.xz*b.wxyzt + a.wxz*b.xyzt - a.xyz*b.wxzt + a.xzt*b.wxyz + a.wxyz*b.xzt - a.wxzt*b.xyz + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
	c.xt =  + a.xt*b.wxyzt + a.wxt*b.xyzt - a.xyt*b.wxzt + a.xzt*b.wxyt + a.wxyt*b.xzt - a.wxzt*b.xyt + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
	c.yz =  + a.yz*b.wxyzt + a.wyz*b.xyzt - a.xyz*b.wyzt + a.yzt*b.wxyz + a.wxyz*b.yzt - a.wyzt*b.xyz + a.xyzt*b.wyz + a.wxyzt*b.yz ; 
	c.yt =  + a.yt*b.wxyzt + a.wyt*b.xyzt - a.xyt*b.wyzt + a.yzt*b.wxyt + a.wxyt*b.yzt - a.wyzt*b.xyt + a.xyzt*b.wyt + a.wxyzt*b.yt ; 
	c.zt =  + a.zt*b.wxyzt + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 

	c.wxy =  + a.wxy*b.wxyzt + a.wxyz*b.wxyt - a.wxyt*b.wxyz + a.wxyzt*b.wxy ; 
	c.wxz =  + a.wxz*b.wxyzt + a.wxyz*b.wxzt - a.wxzt*b.wxyz + a.wxyzt*b.wxz ; 
	c.wxt =  + a.wxt*b.wxyzt + a.wxyt*b.wxzt - a.wxzt*b.wxyt + a.wxyzt*b.wxt ; 
	c.wyz =  + a.wyz*b.wxyzt + a.wxyz*b.wyzt - a.wyzt*b.wxyz + a.wxyzt*b.wyz ; 
	c.wyt =  + a.wyt*b.wxyzt + a.wxyt*b.wyzt - a.wyzt*b.wxyt + a.wxyzt*b.wyt ; 
	c.wzt =  + a.wzt*b.wxyzt + a.wxzt*b.wyzt - a.wyzt*b.wxzt + a.wxyzt*b.wzt ; 
	c.xyz =  + a.xyz*b.wxyzt + a.wxyz*b.xyzt - a.xyzt*b.wxyz + a.wxyzt*b.xyz ; 
	c.xyt =  + a.xyt*b.wxyzt + a.wxyt*b.xyzt - a.xyzt*b.wxyt + a.wxyzt*b.xyt ; 
	c.xzt =  + a.xzt*b.wxyzt + a.wxzt*b.xyzt - a.xyzt*b.wxzt + a.wxyzt*b.xzt ; 
	c.yzt =  + a.yzt*b.wxyzt + a.wyzt*b.xyzt - a.xyzt*b.wyzt + a.wxyzt*b.yzt ; 

	c.wxyz =  + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz ; 
	c.wxyt =  + a.wxyt*b.wxyzt + a.wxyzt*b.wxyt ; 
	c.wxzt =  + a.wxzt*b.wxyzt + a.wxyzt*b.wxzt ; 
	c.wyzt =  + a.wyzt*b.wxyzt + a.wxyzt*b.wyzt ; 
	c.xyzt =  + a.xyzt*b.wxyzt + a.wxyzt*b.xyzt ; 

	c.wxyzt =  + a.wxyzt*b.wxyzt ; 

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Regressive(GA5_4_1 a, GA5_4_1 b) 
{
	GA5_4_1 c;

//	GA5_4_1 I, I_inv;
//	I.wxyzt = 1;
//	I_inv.wxyzt = -1;
//	c = ((a*I_inv)^(b*I_inv))*I;


	c.q =  + a.q*b.wxyzt + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz 
		+ a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt 
		+ a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt
		 + a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz
		 - a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x + a.xyzt*b.w + a.wxyzt*b.q ; 

	c.w =  + a.w*b.wxyzt - a.wx*b.wyzt + a.wy*b.wxzt - a.wz*b.wxyt + a.wt*b.wxyz
		 + a.wxy*b.wzt - a.wxz*b.wyt + a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz
		 + a.wzt*b.wxy - a.wxyz*b.wt + a.wxyt*b.wz - a.wxzt*b.wy + a.wyzt*b.wx + a.wxyzt*b.w ; 
	c.x =  + a.x*b.wxyzt - a.wx*b.xyzt + a.xy*b.wxzt - a.xz*b.wxyt + a.xt*b.wxyz
		 + a.wxy*b.xzt - a.wxz*b.xyt + a.wxt*b.xyz + a.xyz*b.wxt - a.xyt*b.wxz
		 + a.xzt*b.wxy - a.wxyz*b.xt + a.wxyt*b.xz - a.wxzt*b.xy + a.xyzt*b.wx + a.wxyzt*b.x ; 
	c.y =  + a.y*b.wxyzt - a.wy*b.xyzt + a.xy*b.wyzt - a.yz*b.wxyt + a.yt*b.wxyz
		 + a.wxy*b.yzt - a.wyz*b.xyt + a.wyt*b.xyz + a.xyz*b.wyt - a.xyt*b.wyz
		 + a.yzt*b.wxy - a.wxyz*b.yt + a.wxyt*b.yz - a.wyzt*b.xy + a.xyzt*b.wy + a.wxyzt*b.y ; 
	c.z =  + a.z*b.wxyzt - a.wz*b.xyzt + a.xz*b.wyzt - a.yz*b.wxzt + a.zt*b.wxyz
		 + a.wxz*b.yzt - a.wyz*b.xzt + a.wzt*b.xyz + a.xyz*b.wzt - a.xzt*b.wyz
		 + a.yzt*b.wxz - a.wxyz*b.zt + a.wxzt*b.yz - a.wyzt*b.xz + a.xyzt*b.wz + a.wxyzt*b.z ; 
	c.t =  + a.t*b.wxyzt - a.wt*b.xyzt + a.xt*b.wyzt - a.yt*b.wxzt + a.zt*b.wxyt
		 + a.wxt*b.yzt - a.wyt*b.xzt + a.wzt*b.xyt + a.xyt*b.wzt - a.xzt*b.wyt
		 + a.yzt*b.wxt - a.wxyt*b.zt + a.wxzt*b.yt - a.wyzt*b.xt + a.xyzt*b.wt + a.wxyzt*b.t ; 

	c.wx =  + a.wx*b.wxyzt + a.wxy*b.wxzt - a.wxz*b.wxyt + a.wxt*b.wxyz + a.wxyz*b.wxt - a.wxyt*b.wxz + a.wxzt*b.wxy + a.wxyzt*b.wx ; 
	c.wy =  + a.wy*b.wxyzt + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz + a.wxyz*b.wyt - a.wxyt*b.wyz + a.wyzt*b.wxy + a.wxyzt*b.wy ; 
	c.wz =  + a.wz*b.wxyzt + a.wxz*b.wyzt - a.wyz*b.wxzt + a.wzt*b.wxyz + a.wxyz*b.wzt - a.wxzt*b.wyz + a.wyzt*b.wxz + a.wxyzt*b.wz ; 
	c.wt =  + a.wt*b.wxyzt + a.wxt*b.wyzt - a.wyt*b.wxzt + a.wzt*b.wxyt + a.wxyt*b.wzt - a.wxzt*b.wyt + a.wyzt*b.wxt + a.wxyzt*b.wt ; 
	c.xy =  + a.xy*b.wxyzt + a.wxy*b.xyzt - a.xyz*b.wxyt + a.xyt*b.wxyz + a.wxyz*b.xyt - a.wxyt*b.xyz + a.xyzt*b.wxy + a.wxyzt*b.xy ; 
	c.xz =  + a.xz*b.wxyzt + a.wxz*b.xyzt - a.xyz*b.wxzt + a.xzt*b.wxyz + a.wxyz*b.xzt - a.wxzt*b.xyz + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
	c.xt =  + a.xt*b.wxyzt + a.wxt*b.xyzt - a.xyt*b.wxzt + a.xzt*b.wxyt + a.wxyt*b.xzt - a.wxzt*b.xyt + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
	c.yz =  + a.yz*b.wxyzt + a.wyz*b.xyzt - a.xyz*b.wyzt + a.yzt*b.wxyz + a.wxyz*b.yzt - a.wyzt*b.xyz + a.xyzt*b.wyz + a.wxyzt*b.yz ; 
	c.yt =  + a.yt*b.wxyzt + a.wyt*b.xyzt - a.xyt*b.wyzt + a.yzt*b.wxyt + a.wxyt*b.yzt - a.wyzt*b.xyt + a.xyzt*b.wyt + a.wxyzt*b.yt ; 
	c.zt =  + a.zt*b.wxyzt + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 

	c.wxy =  + a.wxy*b.wxyzt - a.wxyz*b.wxyt + a.wxyt*b.wxyz + a.wxyzt*b.wxy ; 
	c.wxz =  + a.wxz*b.wxyzt - a.wxyz*b.wxzt + a.wxzt*b.wxyz + a.wxyzt*b.wxz ; 
	c.wxt =  + a.wxt*b.wxyzt - a.wxyt*b.wxzt + a.wxzt*b.wxyt + a.wxyzt*b.wxt ; 
	c.wyz =  + a.wyz*b.wxyzt - a.wxyz*b.wyzt + a.wyzt*b.wxyz + a.wxyzt*b.wyz ; 
	c.wyt =  + a.wyt*b.wxyzt - a.wxyt*b.wyzt + a.wyzt*b.wxyt + a.wxyzt*b.wyt ; 
	c.wzt =  + a.wzt*b.wxyzt - a.wxzt*b.wyzt + a.wyzt*b.wxzt + a.wxyzt*b.wzt ; 
	c.xyz =  + a.xyz*b.wxyzt - a.wxyz*b.xyzt + a.xyzt*b.wxyz + a.wxyzt*b.xyz ; 
	c.xyt =  + a.xyt*b.wxyzt - a.wxyt*b.xyzt + a.xyzt*b.wxyt + a.wxyzt*b.xyt ; 
	c.xzt =  + a.xzt*b.wxyzt - a.wxzt*b.xyzt + a.xyzt*b.wxzt + a.wxyzt*b.xzt ; 
	c.yzt =  + a.yzt*b.wxyzt - a.wyzt*b.xyzt + a.xyzt*b.wyzt + a.wxyzt*b.yzt ; 

	c.wxyz =  + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz ; 
	c.wxyt =  + a.wxyt*b.wxyzt + a.wxyzt*b.wxyt ; 
	c.wxzt =  + a.wxzt*b.wxyzt + a.wxyzt*b.wxzt ; 
	c.wyzt =  + a.wyzt*b.wxyzt + a.wxyzt*b.wyzt ; 
	c.xyzt =  + a.xyzt*b.wxyzt + a.wxyzt*b.xyzt ; 

	c.wxyzt =  + a.wxyzt*b.wxyzt ; 

	return c;
}


//////////////////////////////////////////////////////

GA5_4_1 RegressiveViaFormula(GA5_4_1 a, GA5_4_1 b) 
{
	GA5_4_1 c;

	GA5_4_1 I, I_inv;
	I.wxyzt = 1;
	I_inv.wxyzt = -1;
	c = Product(Wedge(Product(a,I_inv),Product(b,I_inv)),I);

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 LowerRightViaFormula(GA5_4_1 a, GA5_4_1 b)  // Like a mirror of Wedge along rising diagonal
{

// -Wedge(Blade[i]*xyzt,xyzt*Blade[j])

	GA5_4_1 c;

	GA5_4_1 I,I_inv;
	
	I.wxyzt = 1;	I_inv.wxyzt = -1;

	c = Wedge(Product(a,I_inv),Product(I,b));
	

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Expander(GA5_4_1 a, GA5_4_1 b)
{

	GA5_4_1 c;

//	Expander equation set
	c.q = 0 ; 

	c.w = 0 ; 
	c.x = 0 ; 
	c.y = 0 ; 
	c.z = 0 ; 
	c.t = 0 ; 

	c.wx =  + a.w*b.x - a.x*b.w ; 
	c.wy =  + a.w*b.y - a.y*b.w ; 
	c.wz =  + a.w*b.z - a.z*b.w ; 
	c.wt =  + a.w*b.t - a.t*b.w ; 
	c.xy =  + a.x*b.y - a.y*b.x ; 
	c.xz =  + a.x*b.z - a.z*b.x ; 
	c.xt =  + a.x*b.t - a.t*b.x ; 
	c.yz =  + a.y*b.z - a.z*b.y ; 
	c.yt =  + a.y*b.t - a.t*b.y ; 
	c.zt =  + a.z*b.t - a.t*b.z ; 

	c.wxy =  + a.w*b.xy - a.x*b.wy + a.y*b.wx + a.wx*b.y - a.wy*b.x + a.xy*b.w ; 
	c.wxz =  + a.w*b.xz - a.x*b.wz + a.z*b.wx + a.wx*b.z - a.wz*b.x + a.xz*b.w ; 
	c.wxt =  + a.w*b.xt - a.x*b.wt + a.t*b.wx + a.wx*b.t - a.wt*b.x + a.xt*b.w ; 
	c.wyz =  + a.w*b.yz - a.y*b.wz + a.z*b.wy + a.wy*b.z - a.wz*b.y + a.yz*b.w ; 
	c.wyt =  + a.w*b.yt - a.y*b.wt + a.t*b.wy + a.wy*b.t - a.wt*b.y + a.yt*b.w ; 
	c.wzt =  + a.w*b.zt - a.z*b.wt + a.t*b.wz + a.wz*b.t - a.wt*b.z + a.zt*b.w ; 
	c.xyz =  + a.x*b.yz - a.y*b.xz + a.z*b.xy + a.xy*b.z - a.xz*b.y + a.yz*b.x ; 
	c.xyt =  + a.x*b.yt - a.y*b.xt + a.t*b.xy + a.xy*b.t - a.xt*b.y + a.yt*b.x ; 
	c.xzt =  + a.x*b.zt - a.z*b.xt + a.t*b.xz + a.xz*b.t - a.xt*b.z + a.zt*b.x ; 
	c.yzt =  + a.y*b.zt - a.z*b.yt + a.t*b.yz + a.yz*b.t - a.yt*b.z + a.zt*b.y ; 

	c.wxyz =  + a.w*b.xyz - a.x*b.wyz + a.y*b.wxz - a.z*b.wxy + a.wx*b.yz - a.wy*b.xz 
		+ a.wz*b.xy + a.xy*b.wz - a.xz*b.wy + a.yz*b.wx + a.wxy*b.z - a.wxz*b.y 
		- a.wxt*b.yzt + a.wyz*b.x + a.wyt*b.xzt - a.wzt*b.xyt - a.xyz*b.w 
		- a.xyt*b.wzt + a.xzt*b.wyt - a.yzt*b.wxt ; 
	c.wxyt =  + a.w*b.xyt - a.x*b.wyt + a.y*b.wxt - a.t*b.wxy + a.wx*b.yt - a.wy*b.xt 
		+ a.wt*b.xy + a.xy*b.wt - a.xt*b.wy + a.yt*b.wx + a.wxy*b.t - a.wxz*b.yzt 
		- a.wxt*b.y + a.wyz*b.xzt + a.wyt*b.x - a.wzt*b.xyz - a.xyz*b.wzt 
		- a.xyt*b.w + a.xzt*b.wyz - a.yzt*b.wxz ; 
	c.wxzt =  + a.w*b.xzt - a.x*b.wzt + a.z*b.wxt - a.t*b.wxz + a.wx*b.zt - a.wz*b.xt 
		+ a.wt*b.xz + a.xz*b.wt - a.xt*b.wz + a.zt*b.wx + a.wxy*b.yzt + a.wxz*b.t 
		- a.wxt*b.z - a.wyz*b.xyt + a.wyt*b.xyz + a.wzt*b.x + a.xyz*b.wyt 
		- a.xyt*b.wyz - a.xzt*b.w + a.yzt*b.wxy ; 
	c.wyzt =  + a.w*b.yzt - a.y*b.wzt + a.z*b.wyt - a.t*b.wyz + a.wy*b.zt - a.wz*b.yt 
		+ a.wt*b.yz + a.yz*b.wt - a.yt*b.wz + a.zt*b.wy - a.wxy*b.xzt + a.wxz*b.xyt 
		- a.wxt*b.xyz + a.wyz*b.t - a.wyt*b.z + a.wzt*b.y - a.xyz*b.wxt + a.xyt*b.wxz 
		- a.xzt*b.wxy - a.yzt*b.w ; 
	c.xyzt =  + a.x*b.yzt - a.y*b.xzt + a.z*b.xyt - a.t*b.xyz + a.xy*b.zt - a.xz*b.yt 
		+ a.xt*b.yz + a.yz*b.xt - a.yt*b.xz + a.zt*b.xy + a.wxy*b.wzt - a.wxz*b.wyt 
		+ a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz + a.wzt*b.wxy + a.xyz*b.t 
		- a.xyt*b.z + a.xzt*b.y - a.yzt*b.x ; 

	c.wxyzt =  + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz 
		+ a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt 
		+ a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt 
		+ a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz 
		- a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x + a.xyzt*b.w ; 

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Conserver(GA5_4_1 a, GA5_4_1 b)
{

	GA5_4_1 c;

// Conserver equation set

	c.q =  + a.q*b.q ; 

	c.w =  + a.q*b.w + a.w*b.q ; 
	c.x =  + a.q*b.x + a.x*b.q ; 
	c.y =  + a.q*b.y + a.y*b.q ; 
	c.z =  + a.q*b.z + a.z*b.q ; 
	c.t =  + a.q*b.t + a.t*b.q ; 

	c.wx =  + a.q*b.wx + a.wx*b.q - a.wy*b.xy - a.wz*b.xz + a.wt*b.xt + a.xy*b.wy + a.xz*b.wz - a.xt*b.wt ; 
	c.wy =  + a.q*b.wy + a.wx*b.xy + a.wy*b.q - a.wz*b.yz + a.wt*b.yt - a.xy*b.wx + a.yz*b.wz - a.yt*b.wt ; 
	c.wz =  + a.q*b.wz + a.wx*b.xz + a.wy*b.yz + a.wz*b.q + a.wt*b.zt - a.xz*b.wx - a.yz*b.wy - a.zt*b.wt ; 
	c.wt =  + a.q*b.wt + a.wx*b.xt + a.wy*b.yt + a.wz*b.zt + a.wt*b.q - a.xt*b.wx - a.yt*b.wy - a.zt*b.wz ; 
	c.xy =  + a.q*b.xy - a.wx*b.wy + a.wy*b.wx + a.xy*b.q - a.xz*b.yz + a.xt*b.yt + a.yz*b.xz - a.yt*b.xt ; 
	c.xz =  + a.q*b.xz - a.wx*b.wz + a.wz*b.wx + a.xy*b.yz + a.xz*b.q + a.xt*b.zt - a.yz*b.xy - a.zt*b.xt ; 
	c.xt =  + a.q*b.xt - a.wx*b.wt + a.wt*b.wx + a.xy*b.yt + a.xz*b.zt + a.xt*b.q - a.yt*b.xy - a.zt*b.xz ; 
	c.yz =  + a.q*b.yz - a.wy*b.wz + a.wz*b.wy - a.xy*b.xz + a.xz*b.xy + a.yz*b.q + a.yt*b.zt - a.zt*b.yt ; 
	c.yt =  + a.q*b.yt - a.wy*b.wt + a.wt*b.wy - a.xy*b.xt + a.xt*b.xy + a.yz*b.zt + a.yt*b.q - a.zt*b.yz ; 
	c.zt =  + a.q*b.zt - a.wz*b.wt + a.wt*b.wz - a.xz*b.xt + a.xt*b.xz - a.yz*b.yt + a.yt*b.yz + a.zt*b.q ; 

	c.wxy =  + a.q*b.wxy + a.wz*b.xyz - a.wt*b.xyt - a.xz*b.wyz + a.xt*b.wyt + a.yz*b.wxz - a.yt*b.wxt
		 + a.wxy*b.q - a.wxz*b.yz + a.wxt*b.yt + a.wyz*b.xz - a.wyt*b.xt - a.xyz*b.wz + a.xyt*b.wt ; 
	c.wxz =  + a.q*b.wxz - a.wy*b.xyz - a.wt*b.xzt + a.xy*b.wyz + a.xt*b.wzt - a.yz*b.wxy - a.zt*b.wxt
		 + a.wxy*b.yz + a.wxz*b.q + a.wxt*b.zt - a.wyz*b.xy - a.wzt*b.xt + a.xyz*b.wy + a.xzt*b.wt ; 
	c.wxt =  + a.q*b.wxt - a.wy*b.xyt - a.wz*b.xzt + a.xy*b.wyt + a.xz*b.wzt - a.yt*b.wxy - a.zt*b.wxz
		 + a.wxy*b.yt + a.wxz*b.zt + a.wxt*b.q - a.wyt*b.xy - a.wzt*b.xz + a.xyt*b.wy + a.xzt*b.wz ; 
	c.wyz =  + a.q*b.wyz + a.wx*b.xyz - a.wt*b.yzt - a.xy*b.wxz + a.xz*b.wxy + a.yt*b.wzt - a.zt*b.wyt
		 - a.wxy*b.xz + a.wxz*b.xy + a.wyz*b.q + a.wyt*b.zt - a.wzt*b.yt - a.xyz*b.wx + a.yzt*b.wt ; 
	c.wyt =  + a.q*b.wyt + a.wx*b.xyt - a.wz*b.yzt - a.xy*b.wxt + a.xt*b.wxy + a.yz*b.wzt - a.zt*b.wyz
		 - a.wxy*b.xt + a.wxt*b.xy + a.wyz*b.zt + a.wyt*b.q - a.wzt*b.yz - a.xyt*b.wx + a.yzt*b.wz ; 
	c.wzt =  + a.q*b.wzt + a.wx*b.xzt + a.wy*b.yzt - a.xz*b.wxt + a.xt*b.wxz - a.yz*b.wyt + a.yt*b.wyz
		 - a.wxz*b.xt + a.wxt*b.xz - a.wyz*b.yt + a.wyt*b.yz + a.wzt*b.q - a.xzt*b.wx - a.yzt*b.wy ; 
	c.xyz =  + a.q*b.xyz - a.wx*b.wyz + a.wy*b.wxz - a.wz*b.wxy - a.xt*b.yzt + a.yt*b.xzt - a.zt*b.xyt
		 + a.wxy*b.wz - a.wxz*b.wy + a.wyz*b.wx + a.xyz*b.q + a.xyt*b.zt - a.xzt*b.yt + a.yzt*b.xt ; 
	c.xyt =  + a.q*b.xyt - a.wx*b.wyt + a.wy*b.wxt - a.wt*b.wxy - a.xz*b.yzt + a.yz*b.xzt - a.zt*b.xyz
		 + a.wxy*b.wt - a.wxt*b.wy + a.wyt*b.wx + a.xyz*b.zt + a.xyt*b.q - a.xzt*b.yz + a.yzt*b.xz ; 
	c.xzt =  + a.q*b.xzt - a.wx*b.wzt + a.wz*b.wxt - a.wt*b.wxz + a.xy*b.yzt - a.yz*b.xyt + a.yt*b.xyz
		 + a.wxz*b.wt - a.wxt*b.wz + a.wzt*b.wx - a.xyz*b.yt + a.xyt*b.yz + a.xzt*b.q - a.yzt*b.xy ; 
	c.yzt =  + a.q*b.yzt - a.wy*b.wzt + a.wz*b.wyt - a.wt*b.wyz - a.xy*b.xzt + a.xz*b.xyt - a.xt*b.xyz
		 + a.wyz*b.wt - a.wyt*b.wz + a.wzt*b.wy + a.xyz*b.xt - a.xyt*b.xz + a.xzt*b.xy + a.yzt*b.q ; 

	c.wxyz =  + a.q*b.wxyz + a.wt*b.xyzt - a.xt*b.wyzt + a.yt*b.wxzt - a.zt*b.wxyt
		 + a.wxyz*b.q + a.wxyt*b.zt - a.wxzt*b.yt + a.wyzt*b.xt - a.xyzt*b.wt ; 
	c.wxyt =  + a.q*b.wxyt + a.wz*b.xyzt - a.xz*b.wyzt + a.yz*b.wxzt - a.zt*b.wxyz
		 + a.wxyz*b.zt + a.wxyt*b.q - a.wxzt*b.yz + a.wyzt*b.xz - a.xyzt*b.wz ; 
	c.wxzt =  + a.q*b.wxzt - a.wy*b.xyzt + a.xy*b.wyzt - a.yz*b.wxyt + a.yt*b.wxyz
		 - a.wxyz*b.yt + a.wxyt*b.yz + a.wxzt*b.q - a.wyzt*b.xy + a.xyzt*b.wy ; 
	c.wyzt =  + a.q*b.wyzt + a.wx*b.xyzt - a.xy*b.wxzt + a.xz*b.wxyt - a.xt*b.wxyz
		 + a.wxyz*b.xt - a.wxyt*b.xz + a.wxzt*b.xy + a.wyzt*b.q - a.xyzt*b.wx ; 
	c.xyzt =  + a.q*b.xyzt - a.wx*b.wyzt + a.wy*b.wxzt - a.wz*b.wxyt + a.wt*b.wxyz
		 - a.wxyz*b.wt + a.wxyt*b.wz - a.wxzt*b.wy + a.wyzt*b.wx + a.xyzt*b.q ; 

	c.wxyzt =  + a.q*b.wxyzt + a.wxyzt*b.q ; 

	return c;
}


//////////////////////////////////////////////////////

GA5_4_1 Shrinker(GA5_4_1 a, GA5_4_1 b)
{

	GA5_4_1 c;

// Shrinker equation set

	c.q =  + a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z - a.t*b.t - a.wx*b.wx - a.wy*b.wy
		 - a.wz*b.wz + a.wt*b.wt - a.xy*b.xy - a.xz*b.xz + a.xt*b.xt - a.yz*b.yz + a.yt*b.yt + a.zt*b.zt
		 - a.wxy*b.wxy - a.wxz*b.wxz + a.wxt*b.wxt - a.wyz*b.wyz + a.wyt*b.wyt + a.wzt*b.wzt - a.xyz*b.xyz
		 + a.xyt*b.xyt + a.xzt*b.xzt + a.yzt*b.yzt + a.wxyz*b.wxyz - a.wxyt*b.wxyt - a.wxzt*b.wxzt
		 - a.wyzt*b.wyzt - a.xyzt*b.xyzt - a.wxyzt*b.wxyzt ; 

	c.w =  - a.x*b.wx - a.y*b.wy - a.z*b.wz + a.t*b.wt + a.wx*b.x + a.wy*b.y + a.wz*b.z - a.wt*b.t
		 - a.xy*b.wxy - a.xz*b.wxz + a.xt*b.wxt - a.yz*b.wyz + a.yt*b.wyt + a.zt*b.wzt
		 - a.wxy*b.xy - a.wxz*b.xz + a.wxt*b.xt - a.wyz*b.yz + a.wyt*b.yt + a.wzt*b.zt
		 + a.xyz*b.wxyz - a.xyt*b.wxyt - a.xzt*b.wxzt - a.yzt*b.wyzt - a.wxyz*b.xyz
		 + a.wxyt*b.xyt + a.wxzt*b.xzt + a.wyzt*b.yzt - a.xyzt*b.wxyzt - a.wxyzt*b.xyzt ; 
	c.x =  + a.w*b.wx - a.y*b.xy - a.z*b.xz + a.t*b.xt - a.wx*b.w + a.wy*b.wxy + a.wz*b.wxz
		 - a.wt*b.wxt + a.xy*b.y + a.xz*b.z - a.xt*b.t - a.yz*b.xyz + a.yt*b.xyt + a.zt*b.xzt
		 + a.wxy*b.wy + a.wxz*b.wz - a.wxt*b.wt - a.wyz*b.wxyz + a.wyt*b.wxyt + a.wzt*b.wxzt
		 - a.xyz*b.yz + a.xyt*b.yt + a.xzt*b.zt - a.yzt*b.xyzt + a.wxyz*b.wyz - a.wxyt*b.wyt
		 - a.wxzt*b.wzt + a.wyzt*b.wxyzt + a.xyzt*b.yzt + a.wxyzt*b.wyzt ; 
	c.y =  + a.w*b.wy + a.x*b.xy - a.z*b.yz + a.t*b.yt - a.wx*b.wxy - a.wy*b.w + a.wz*b.wyz
		 - a.wt*b.wyt - a.xy*b.x + a.xz*b.xyz - a.xt*b.xyt + a.yz*b.z - a.yt*b.t + a.zt*b.yzt
		 - a.wxy*b.wx + a.wxz*b.wxyz - a.wxt*b.wxyt + a.wyz*b.wz - a.wyt*b.wt + a.wzt*b.wyzt
		 + a.xyz*b.xz - a.xyt*b.xt + a.xzt*b.xyzt + a.yzt*b.zt - a.wxyz*b.wxz + a.wxyt*b.wxt
		 - a.wxzt*b.wxyzt - a.wyzt*b.wzt - a.xyzt*b.xzt - a.wxyzt*b.wxzt ; 
	c.z =  + a.w*b.wz + a.x*b.xz + a.y*b.yz + a.t*b.zt - a.wx*b.wxz - a.wy*b.wyz - a.wz*b.w
		 - a.wt*b.wzt - a.xy*b.xyz - a.xz*b.x - a.xt*b.xzt - a.yz*b.y - a.yt*b.yzt - a.zt*b.t
		 - a.wxy*b.wxyz - a.wxz*b.wx - a.wxt*b.wxzt - a.wyz*b.wy - a.wyt*b.wyzt - a.wzt*b.wt
		 - a.xyz*b.xy - a.xyt*b.xyzt - a.xzt*b.xt - a.yzt*b.yt + a.wxyz*b.wxy + a.wxyt*b.wxyzt
		 + a.wxzt*b.wxt + a.wyzt*b.wyt + a.xyzt*b.xyt + a.wxyzt*b.wxyt ; 
	c.t =  + a.w*b.wt + a.x*b.xt + a.y*b.yt + a.z*b.zt - a.wx*b.wxt - a.wy*b.wyt - a.wz*b.wzt
		 - a.wt*b.w - a.xy*b.xyt - a.xz*b.xzt - a.xt*b.x - a.yz*b.yzt - a.yt*b.y - a.zt*b.z
		 - a.wxy*b.wxyt - a.wxz*b.wxzt - a.wxt*b.wx - a.wyz*b.wyzt - a.wyt*b.wy - a.wzt*b.wz
		 - a.xyz*b.xyzt - a.xyt*b.xy - a.xzt*b.xz - a.yzt*b.yz + a.wxyz*b.wxyzt + a.wxyt*b.wxy
		 + a.wxzt*b.wxz + a.wyzt*b.wyz + a.xyzt*b.xyz + a.wxyzt*b.wxyz ; 

	c.wx =  + a.y*b.wxy + a.z*b.wxz - a.t*b.wxt - a.yz*b.wxyz + a.yt*b.wxyt + a.zt*b.wxzt + a.wxy*b.y
		 + a.wxz*b.z - a.wxt*b.t - a.wyz*b.xyz + a.wyt*b.xyt + a.wzt*b.xzt + a.xyz*b.wyz - a.xyt*b.wyt
		 - a.xzt*b.wzt + a.yzt*b.wxyzt - a.wxyz*b.yz + a.wxyt*b.yt + a.wxzt*b.zt - a.wyzt*b.xyzt
		 + a.xyzt*b.wyzt + a.wxyzt*b.yzt ; 
	c.wy =  - a.x*b.wxy + a.z*b.wyz - a.t*b.wyt + a.xz*b.wxyz - a.xt*b.wxyt + a.zt*b.wyzt - a.wxy*b.x
		 + a.wxz*b.xyz - a.wxt*b.xyt + a.wyz*b.z - a.wyt*b.t + a.wzt*b.yzt - a.xyz*b.wxz + a.xyt*b.wxt
		 - a.xzt*b.wxyzt - a.yzt*b.wzt + a.wxyz*b.xz - a.wxyt*b.xt + a.wxzt*b.xyzt + a.wyzt*b.zt
		 - a.xyzt*b.wxzt - a.wxyzt*b.xzt ; 
	c.wz =  - a.x*b.wxz - a.y*b.wyz - a.t*b.wzt - a.xy*b.wxyz - a.xt*b.wxzt - a.yt*b.wyzt - a.wxy*b.xyz
		 - a.wxz*b.x - a.wxt*b.xzt - a.wyz*b.y - a.wyt*b.yzt - a.wzt*b.t + a.xyz*b.wxy + a.xyt*b.wxyzt
		 + a.xzt*b.wxt + a.yzt*b.wyt - a.wxyz*b.xy - a.wxyt*b.xyzt - a.wxzt*b.xt - a.wyzt*b.yt
		 + a.xyzt*b.wxyt + a.wxyzt*b.xyt ; 
	c.wt =  - a.x*b.wxt - a.y*b.wyt - a.z*b.wzt - a.xy*b.wxyt - a.xz*b.wxzt - a.yz*b.wyzt - a.wxy*b.xyt
		 - a.wxz*b.xzt - a.wxt*b.x - a.wyz*b.yzt - a.wyt*b.y - a.wzt*b.z + a.xyz*b.wxyzt + a.xyt*b.wxy
		 + a.xzt*b.wxz + a.yzt*b.wyz - a.wxyz*b.xyzt - a.wxyt*b.xy - a.wxzt*b.xz - a.wyzt*b.yz
		 + a.xyzt*b.wxyz + a.wxyzt*b.xyz ; 
	c.xy =  + a.w*b.wxy + a.z*b.xyz - a.t*b.xyt - a.wz*b.wxyz + a.wt*b.wxyt + a.zt*b.xyzt + a.wxy*b.w
		 - a.wxz*b.wyz + a.wxt*b.wyt + a.wyz*b.wxz - a.wyt*b.wxt + a.wzt*b.wxyzt + a.xyz*b.z - a.xyt*b.t
		 + a.xzt*b.yzt - a.yzt*b.xzt - a.wxyz*b.wz + a.wxyt*b.wt - a.wxzt*b.wyzt + a.wyzt*b.wxzt
		 + a.xyzt*b.zt + a.wxyzt*b.wzt ; 
	c.xz =  + a.w*b.wxz - a.y*b.xyz - a.t*b.xzt + a.wy*b.wxyz + a.wt*b.wxzt - a.yt*b.xyzt + a.wxy*b.wyz
		 + a.wxz*b.w + a.wxt*b.wzt - a.wyz*b.wxy - a.wyt*b.wxyzt - a.wzt*b.wxt - a.xyz*b.y - a.xyt*b.yzt
		 - a.xzt*b.t + a.yzt*b.xyt + a.wxyz*b.wy + a.wxyt*b.wyzt + a.wxzt*b.wt - a.wyzt*b.wxyt
		 - a.xyzt*b.yt - a.wxyzt*b.wyt ; 
	c.xt =  + a.w*b.wxt - a.y*b.xyt - a.z*b.xzt + a.wy*b.wxyt + a.wz*b.wxzt - a.yz*b.xyzt + a.wxy*b.wyt
		 + a.wxz*b.wzt + a.wxt*b.w - a.wyz*b.wxyzt - a.wyt*b.wxy - a.wzt*b.wxz - a.xyz*b.yzt - a.xyt*b.y
		 - a.xzt*b.z + a.yzt*b.xyz + a.wxyz*b.wyzt + a.wxyt*b.wy + a.wxzt*b.wz - a.wyzt*b.wxyz
		 - a.xyzt*b.yz - a.wxyzt*b.wyz ; 
	c.yz =  + a.w*b.wyz + a.x*b.xyz - a.t*b.yzt - a.wx*b.wxyz + a.wt*b.wyzt + a.xt*b.xyzt - a.wxy*b.wxz
		 + a.wxz*b.wxy + a.wxt*b.wxyzt + a.wyz*b.w + a.wyt*b.wzt - a.wzt*b.wyt + a.xyz*b.x + a.xyt*b.xzt
		 - a.xzt*b.xyt - a.yzt*b.t - a.wxyz*b.wx - a.wxyt*b.wxzt + a.wxzt*b.wxyt + a.wyzt*b.wt
		 + a.xyzt*b.xt + a.wxyzt*b.wxt ; 
	c.yt =  + a.w*b.wyt + a.x*b.xyt - a.z*b.yzt - a.wx*b.wxyt + a.wz*b.wyzt + a.xz*b.xyzt - a.wxy*b.wxt
		 + a.wxz*b.wxyzt + a.wxt*b.wxy + a.wyz*b.wzt + a.wyt*b.w - a.wzt*b.wyz + a.xyz*b.xzt + a.xyt*b.x
		 - a.xzt*b.xyz - a.yzt*b.z - a.wxyz*b.wxzt - a.wxyt*b.wx + a.wxzt*b.wxyz + a.wyzt*b.wz
		 + a.xyzt*b.xz + a.wxyzt*b.wxz ; 
	c.zt =  + a.w*b.wzt + a.x*b.xzt + a.y*b.yzt - a.wx*b.wxzt - a.wy*b.wyzt - a.xy*b.xyzt - a.wxy*b.wxyzt
		 - a.wxz*b.wxt + a.wxt*b.wxz - a.wyz*b.wyt + a.wyt*b.wyz + a.wzt*b.w - a.xyz*b.xyt + a.xyt*b.xyz
		 + a.xzt*b.x + a.yzt*b.y + a.wxyz*b.wxyt - a.wxyt*b.wxyz - a.wxzt*b.wx - a.wyzt*b.wy
		 - a.xyzt*b.xy - a.wxyzt*b.wxy ; 

	c.wxy =  - a.z*b.wxyz + a.t*b.wxyt + a.zt*b.wxyzt + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt
		 + a.wxyz*b.z - a.wxyt*b.t + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 
	c.wxz =  + a.y*b.wxyz + a.t*b.wxzt - a.yt*b.wxyzt - a.wyt*b.xyzt + a.xyt*b.wyzt - a.yzt*b.wxyt
		 - a.wxyz*b.y - a.wxyt*b.yzt - a.wxzt*b.t + a.wyzt*b.xyt - a.xyzt*b.wyt - a.wxyzt*b.yt ; 
	c.wxt =  + a.y*b.wxyt + a.z*b.wxzt - a.yz*b.wxyzt - a.wyz*b.xyzt + a.xyz*b.wyzt - a.yzt*b.wxyz
		 - a.wxyz*b.yzt - a.wxyt*b.y - a.wxzt*b.z + a.wyzt*b.xyz - a.xyzt*b.wyz - a.wxyzt*b.yz ; 
	c.wyz =  - a.x*b.wxyz + a.t*b.wyzt + a.xt*b.wxyzt + a.wxt*b.xyzt - a.xyt*b.wxzt + a.xzt*b.wxyt
		 + a.wxyz*b.x + a.wxyt*b.xzt - a.wxzt*b.xyt - a.wyzt*b.t + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
	c.wyt =  - a.x*b.wxyt + a.z*b.wyzt + a.xz*b.wxyzt + a.wxz*b.xyzt - a.xyz*b.wxzt + a.xzt*b.wxyz
		 + a.wxyz*b.xzt + a.wxyt*b.x - a.wxzt*b.xyz - a.wyzt*b.z + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
	c.wzt =  - a.x*b.wxzt - a.y*b.wyzt - a.xy*b.wxyzt - a.wxy*b.xyzt + a.xyz*b.wxyt - a.xyt*b.wxyz
		 - a.wxyz*b.xyt + a.wxyt*b.xyz + a.wxzt*b.x + a.wyzt*b.y - a.xyzt*b.wxy - a.wxyzt*b.xy ; 
	c.xyz =  + a.w*b.wxyz + a.t*b.xyzt - a.wt*b.wxyzt - a.wxt*b.wyzt + a.wyt*b.wxzt - a.wzt*b.wxyt
		 - a.wxyz*b.w - a.wxyt*b.wzt + a.wxzt*b.wyt - a.wyzt*b.wxt - a.xyzt*b.t - a.wxyzt*b.wt ; 
	c.xyt =  + a.w*b.wxyt + a.z*b.xyzt - a.wz*b.wxyzt - a.wxz*b.wyzt + a.wyz*b.wxzt - a.wzt*b.wxyz
		 - a.wxyz*b.wzt - a.wxyt*b.w + a.wxzt*b.wyz - a.wyzt*b.wxz - a.xyzt*b.z - a.wxyzt*b.wz ; 
	c.xzt =  + a.w*b.wxzt - a.y*b.xyzt + a.wy*b.wxyzt + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz
		 + a.wxyz*b.wyt - a.wxyt*b.wyz - a.wxzt*b.w + a.wyzt*b.wxy + a.xyzt*b.y + a.wxyzt*b.wy ; 
	c.yzt =  + a.w*b.wyzt + a.x*b.xyzt - a.wx*b.wxyzt - a.wxy*b.wxzt + a.wxz*b.wxyt - a.wxt*b.wxyz
		 - a.wxyz*b.wxt + a.wxyt*b.wxz - a.wxzt*b.wxy - a.wyzt*b.w - a.xyzt*b.x - a.wxyzt*b.wx ; 

	c.wxyz =  - a.t*b.wxyzt - a.wxyzt*b.t ; 
	c.wxyt =  - a.z*b.wxyzt - a.wxyzt*b.z ; 
	c.wxzt =  + a.y*b.wxyzt + a.wxyzt*b.y ; 
	c.wyzt =  - a.x*b.wxyzt - a.wxyzt*b.x ; 
	c.xyzt =  + a.w*b.wxyzt + a.wxyzt*b.w ; 

	c.wxyzt = 0 ; 

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Symmetric(GA5_4_1 a,  GA5_4_1 b) {

	GA5_4_1 c;

//	c = a*b/2 + b*a/2;


// Symmetric Product
c.q =  + a.q*b.q + a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z - a.t*b.t - a.wx*b.wx - a.wy*b.wy - a.wz*b.wz + a.wt*b.wt - a.xy*b.xy - a.xz*b.xz + a.xt*b.xt - a.yz*b.yz + a.yt*b.yt + a.zt*b.zt - a.wxy*b.wxy - a.wxz*b.wxz + a.wxt*b.wxt - a.wyz*b.wyz + a.wyt*b.wyt + a.wzt*b.wzt - a.xyz*b.xyz + a.xyt*b.xyt + a.xzt*b.xzt + a.yzt*b.yzt + a.wxyz*b.wxyz - a.wxyt*b.wxyt - a.wxzt*b.wxzt - a.wyzt*b.wyzt - a.xyzt*b.xyzt - a.wxyzt*b.wxyzt ; 
c.w =  + a.q*b.w + a.w*b.q - a.xy*b.wxy - a.xz*b.wxz + a.xt*b.wxt - a.yz*b.wyz + a.yt*b.wyt + a.zt*b.wzt - a.wxy*b.xy - a.wxz*b.xz + a.wxt*b.xt - a.wyz*b.yz + a.wyt*b.yt + a.wzt*b.zt - a.xyzt*b.wxyzt - a.wxyzt*b.xyzt ; 
c.x =  + a.q*b.x + a.x*b.q + a.wy*b.wxy + a.wz*b.wxz - a.wt*b.wxt - a.yz*b.xyz + a.yt*b.xyt + a.zt*b.xzt + a.wxy*b.wy + a.wxz*b.wz - a.wxt*b.wt - a.xyz*b.yz + a.xyt*b.yt + a.xzt*b.zt + a.wyzt*b.wxyzt + a.wxyzt*b.wyzt ; 
c.y =  + a.q*b.y + a.y*b.q - a.wx*b.wxy + a.wz*b.wyz - a.wt*b.wyt + a.xz*b.xyz - a.xt*b.xyt + a.zt*b.yzt - a.wxy*b.wx + a.wyz*b.wz - a.wyt*b.wt + a.xyz*b.xz - a.xyt*b.xt + a.yzt*b.zt - a.wxzt*b.wxyzt - a.wxyzt*b.wxzt ; 
c.z =  + a.q*b.z + a.z*b.q - a.wx*b.wxz - a.wy*b.wyz - a.wt*b.wzt - a.xy*b.xyz - a.xt*b.xzt - a.yt*b.yzt - a.wxz*b.wx - a.wyz*b.wy - a.wzt*b.wt - a.xyz*b.xy - a.xzt*b.xt - a.yzt*b.yt + a.wxyt*b.wxyzt + a.wxyzt*b.wxyt ; 
c.t =  + a.q*b.t + a.t*b.q - a.wx*b.wxt - a.wy*b.wyt - a.wz*b.wzt - a.xy*b.xyt - a.xz*b.xzt - a.yz*b.yzt - a.wxt*b.wx - a.wyt*b.wy - a.wzt*b.wz - a.xyt*b.xy - a.xzt*b.xz - a.yzt*b.yz + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz ; 
c.wx =  + a.q*b.wx + a.y*b.wxy + a.z*b.wxz - a.t*b.wxt + a.wx*b.q - a.yz*b.wxyz + a.yt*b.wxyt + a.zt*b.wxzt + a.wxy*b.y + a.wxz*b.z - a.wxt*b.t + a.yzt*b.wxyzt - a.wxyz*b.yz + a.wxyt*b.yt + a.wxzt*b.zt + a.wxyzt*b.yzt ; 
c.wy =  + a.q*b.wy - a.x*b.wxy + a.z*b.wyz - a.t*b.wyt + a.wy*b.q + a.xz*b.wxyz - a.xt*b.wxyt + a.zt*b.wyzt - a.wxy*b.x + a.wyz*b.z - a.wyt*b.t - a.xzt*b.wxyzt + a.wxyz*b.xz - a.wxyt*b.xt + a.wyzt*b.zt - a.wxyzt*b.xzt ; 
c.wz =  + a.q*b.wz - a.x*b.wxz - a.y*b.wyz - a.t*b.wzt + a.wz*b.q - a.xy*b.wxyz - a.xt*b.wxzt - a.yt*b.wyzt - a.wxz*b.x - a.wyz*b.y - a.wzt*b.t + a.xyt*b.wxyzt - a.wxyz*b.xy - a.wxzt*b.xt - a.wyzt*b.yt + a.wxyzt*b.xyt ; 
c.wt =  + a.q*b.wt - a.x*b.wxt - a.y*b.wyt - a.z*b.wzt + a.wt*b.q - a.xy*b.wxyt - a.xz*b.wxzt - a.yz*b.wyzt - a.wxt*b.x - a.wyt*b.y - a.wzt*b.z + a.xyz*b.wxyzt - a.wxyt*b.xy - a.wxzt*b.xz - a.wyzt*b.yz + a.wxyzt*b.xyz ; 
c.xy =  + a.q*b.xy + a.w*b.wxy + a.z*b.xyz - a.t*b.xyt - a.wz*b.wxyz + a.wt*b.wxyt + a.xy*b.q + a.zt*b.xyzt + a.wxy*b.w + a.wzt*b.wxyzt + a.xyz*b.z - a.xyt*b.t - a.wxyz*b.wz + a.wxyt*b.wt + a.xyzt*b.zt + a.wxyzt*b.wzt ; 
c.xz =  + a.q*b.xz + a.w*b.wxz - a.y*b.xyz - a.t*b.xzt + a.wy*b.wxyz + a.wt*b.wxzt + a.xz*b.q - a.yt*b.xyzt + a.wxz*b.w - a.wyt*b.wxyzt - a.xyz*b.y - a.xzt*b.t + a.wxyz*b.wy + a.wxzt*b.wt - a.xyzt*b.yt - a.wxyzt*b.wyt ; 
c.xt =  + a.q*b.xt + a.w*b.wxt - a.y*b.xyt - a.z*b.xzt + a.wy*b.wxyt + a.wz*b.wxzt + a.xt*b.q - a.yz*b.xyzt + a.wxt*b.w - a.wyz*b.wxyzt - a.xyt*b.y - a.xzt*b.z + a.wxyt*b.wy + a.wxzt*b.wz - a.xyzt*b.yz - a.wxyzt*b.wyz ; 
c.yz =  + a.q*b.yz + a.w*b.wyz + a.x*b.xyz - a.t*b.yzt - a.wx*b.wxyz + a.wt*b.wyzt + a.xt*b.xyzt + a.yz*b.q + a.wxt*b.wxyzt + a.wyz*b.w + a.xyz*b.x - a.yzt*b.t - a.wxyz*b.wx + a.wyzt*b.wt + a.xyzt*b.xt + a.wxyzt*b.wxt ; 
c.yt =  + a.q*b.yt + a.w*b.wyt + a.x*b.xyt - a.z*b.yzt - a.wx*b.wxyt + a.wz*b.wyzt + a.xz*b.xyzt + a.yt*b.q + a.wxz*b.wxyzt + a.wyt*b.w + a.xyt*b.x - a.yzt*b.z - a.wxyt*b.wx + a.wyzt*b.wz + a.xyzt*b.xz + a.wxyzt*b.wxz ; 
c.zt =  + a.q*b.zt + a.w*b.wzt + a.x*b.xzt + a.y*b.yzt - a.wx*b.wxzt - a.wy*b.wyzt - a.xy*b.xyzt + a.zt*b.q - a.wxy*b.wxyzt + a.wzt*b.w + a.xzt*b.x + a.yzt*b.y - a.wxzt*b.wx - a.wyzt*b.wy - a.xyzt*b.xy - a.wxyzt*b.wxy ; 
c.wxy =  + a.q*b.wxy + a.w*b.xy - a.x*b.wy + a.y*b.wx + a.wx*b.y - a.wy*b.x + a.xy*b.w + a.zt*b.wxyzt + a.wxy*b.q + a.wzt*b.xyzt - a.xzt*b.wyzt + a.yzt*b.wxzt + a.wxzt*b.yzt - a.wyzt*b.xzt + a.xyzt*b.wzt + a.wxyzt*b.zt ; 
c.wxz =  + a.q*b.wxz + a.w*b.xz - a.x*b.wz + a.z*b.wx + a.wx*b.z - a.wz*b.x + a.xz*b.w - a.yt*b.wxyzt + a.wxz*b.q - a.wyt*b.xyzt + a.xyt*b.wyzt - a.yzt*b.wxyt - a.wxyt*b.yzt + a.wyzt*b.xyt - a.xyzt*b.wyt - a.wxyzt*b.yt ; 
c.wxt =  + a.q*b.wxt + a.w*b.xt - a.x*b.wt + a.t*b.wx + a.wx*b.t - a.wt*b.x + a.xt*b.w - a.yz*b.wxyzt + a.wxt*b.q - a.wyz*b.xyzt + a.xyz*b.wyzt - a.yzt*b.wxyz - a.wxyz*b.yzt + a.wyzt*b.xyz - a.xyzt*b.wyz - a.wxyzt*b.yz ; 
c.wyz =  + a.q*b.wyz + a.w*b.yz - a.y*b.wz + a.z*b.wy + a.wy*b.z - a.wz*b.y + a.xt*b.wxyzt + a.yz*b.w + a.wxt*b.xyzt + a.wyz*b.q - a.xyt*b.wxzt + a.xzt*b.wxyt + a.wxyt*b.xzt - a.wxzt*b.xyt + a.xyzt*b.wxt + a.wxyzt*b.xt ; 
c.wyt =  + a.q*b.wyt + a.w*b.yt - a.y*b.wt + a.t*b.wy + a.wy*b.t - a.wt*b.y + a.xz*b.wxyzt + a.yt*b.w + a.wxz*b.xyzt + a.wyt*b.q - a.xyz*b.wxzt + a.xzt*b.wxyz + a.wxyz*b.xzt - a.wxzt*b.xyz + a.xyzt*b.wxz + a.wxyzt*b.xz ; 
c.wzt =  + a.q*b.wzt + a.w*b.zt - a.z*b.wt + a.t*b.wz + a.wz*b.t - a.wt*b.z - a.xy*b.wxyzt + a.zt*b.w - a.wxy*b.xyzt + a.wzt*b.q + a.xyz*b.wxyt - a.xyt*b.wxyz - a.wxyz*b.xyt + a.wxyt*b.xyz - a.xyzt*b.wxy - a.wxyzt*b.xy ; 
c.xyz =  + a.q*b.xyz + a.x*b.yz - a.y*b.xz + a.z*b.xy - a.wt*b.wxyzt + a.xy*b.z - a.xz*b.y + a.yz*b.x - a.wxt*b.wyzt + a.wyt*b.wxzt - a.wzt*b.wxyt + a.xyz*b.q - a.wxyt*b.wzt + a.wxzt*b.wyt - a.wyzt*b.wxt - a.wxyzt*b.wt ; 
c.xyt =  + a.q*b.xyt + a.x*b.yt - a.y*b.xt + a.t*b.xy - a.wz*b.wxyzt + a.xy*b.t - a.xt*b.y + a.yt*b.x - a.wxz*b.wyzt + a.wyz*b.wxzt - a.wzt*b.wxyz + a.xyt*b.q - a.wxyz*b.wzt + a.wxzt*b.wyz - a.wyzt*b.wxz - a.wxyzt*b.wz ; 
c.xzt =  + a.q*b.xzt + a.x*b.zt - a.z*b.xt + a.t*b.xz + a.wy*b.wxyzt + a.xz*b.t - a.xt*b.z + a.zt*b.x + a.wxy*b.wyzt - a.wyz*b.wxyt + a.wyt*b.wxyz + a.xzt*b.q + a.wxyz*b.wyt - a.wxyt*b.wyz + a.wyzt*b.wxy + a.wxyzt*b.wy ; 
c.yzt =  + a.q*b.yzt + a.y*b.zt - a.z*b.yt + a.t*b.yz - a.wx*b.wxyzt + a.yz*b.t - a.yt*b.z + a.zt*b.y - a.wxy*b.wxzt + a.wxz*b.wxyt - a.wxt*b.wxyz + a.yzt*b.q - a.wxyz*b.wxt + a.wxyt*b.wxz - a.wxzt*b.wxy - a.wxyzt*b.wx ; 
c.wxyz =  + a.q*b.wxyz - a.t*b.wxyzt + a.wx*b.yz - a.wy*b.xz + a.wz*b.xy + a.xy*b.wz - a.xz*b.wy + a.yz*b.wx - a.wxt*b.yzt + a.wyt*b.xzt - a.wzt*b.xyt - a.xyt*b.wzt + a.xzt*b.wyt - a.yzt*b.wxt + a.wxyz*b.q - a.wxyzt*b.t ; 
c.wxyt =  + a.q*b.wxyt - a.z*b.wxyzt + a.wx*b.yt - a.wy*b.xt + a.wt*b.xy + a.xy*b.wt - a.xt*b.wy + a.yt*b.wx - a.wxz*b.yzt + a.wyz*b.xzt - a.wzt*b.xyz - a.xyz*b.wzt + a.xzt*b.wyz - a.yzt*b.wxz + a.wxyt*b.q - a.wxyzt*b.z ; 
c.wxzt =  + a.q*b.wxzt + a.y*b.wxyzt + a.wx*b.zt - a.wz*b.xt + a.wt*b.xz + a.xz*b.wt - a.xt*b.wz + a.zt*b.wx + a.wxy*b.yzt - a.wyz*b.xyt + a.wyt*b.xyz + a.xyz*b.wyt - a.xyt*b.wyz + a.yzt*b.wxy + a.wxzt*b.q + a.wxyzt*b.y ; 
c.wyzt =  + a.q*b.wyzt - a.x*b.wxyzt + a.wy*b.zt - a.wz*b.yt + a.wt*b.yz + a.yz*b.wt - a.yt*b.wz + a.zt*b.wy - a.wxy*b.xzt + a.wxz*b.xyt - a.wxt*b.xyz - a.xyz*b.wxt + a.xyt*b.wxz - a.xzt*b.wxy + a.wyzt*b.q - a.wxyzt*b.x ; 
c.xyzt =  + a.q*b.xyzt + a.w*b.wxyzt + a.xy*b.zt - a.xz*b.yt + a.xt*b.yz + a.yz*b.xt - a.yt*b.xz + a.zt*b.xy + a.wxy*b.wzt - a.wxz*b.wyt + a.wxt*b.wyz + a.wyz*b.wxt - a.wyt*b.wxz + a.wzt*b.wxy + a.xyzt*b.q + a.wxyzt*b.w ; 
c.wxyzt =  + a.q*b.wxyzt + a.w*b.xyzt - a.x*b.wyzt + a.y*b.wxzt - a.z*b.wxyt + a.t*b.wxyz + a.wx*b.yzt - a.wy*b.xzt + a.wz*b.xyt - a.wt*b.xyz + a.xy*b.wzt - a.xz*b.wyt + a.xt*b.wyz + a.yz*b.wxt - a.yt*b.wxz + a.zt*b.wxy + a.wxy*b.zt - a.wxz*b.yt + a.wxt*b.yz + a.wyz*b.xt - a.wyt*b.xz + a.wzt*b.xy - a.xyz*b.wt + a.xyt*b.wz - a.xzt*b.wy + a.yzt*b.wx + a.wxyz*b.t - a.wxyt*b.z + a.wxzt*b.y - a.wyzt*b.x + a.xyzt*b.w + a.wxyzt*b.q ; 

	return c;
}



//////////////////////////////////////////////////////

GA5_4_1 AntiSymmetric(GA5_4_1 a, GA5_4_1 b) {

	GA5_4_1 c;

//	c = a*b/2 - b*a/2;


// AntiSymmetric Product
c.q =  0; 
c.w =  - a.x*b.wx - a.y*b.wy - a.z*b.wz + a.t*b.wt + a.wx*b.x + a.wy*b.y + a.wz*b.z - a.wt*b.t + a.xyz*b.wxyz - a.xyt*b.wxyt - a.xzt*b.wxzt - a.yzt*b.wyzt - a.wxyz*b.xyz + a.wxyt*b.xyt + a.wxzt*b.xzt + a.wyzt*b.yzt ; 
c.x =  + a.w*b.wx - a.y*b.xy - a.z*b.xz + a.t*b.xt - a.wx*b.w + a.xy*b.y + a.xz*b.z - a.xt*b.t - a.wyz*b.wxyz + a.wyt*b.wxyt + a.wzt*b.wxzt - a.yzt*b.xyzt + a.wxyz*b.wyz - a.wxyt*b.wyt - a.wxzt*b.wzt + a.xyzt*b.yzt ; 
c.y =  + a.w*b.wy + a.x*b.xy - a.z*b.yz + a.t*b.yt - a.wy*b.w - a.xy*b.x + a.yz*b.z - a.yt*b.t + a.wxz*b.wxyz - a.wxt*b.wxyt + a.wzt*b.wyzt + a.xzt*b.xyzt - a.wxyz*b.wxz + a.wxyt*b.wxt - a.wyzt*b.wzt - a.xyzt*b.xzt ; 
c.z =  + a.w*b.wz + a.x*b.xz + a.y*b.yz + a.t*b.zt - a.wz*b.w - a.xz*b.x - a.yz*b.y - a.zt*b.t - a.wxy*b.wxyz - a.wxt*b.wxzt - a.wyt*b.wyzt - a.xyt*b.xyzt + a.wxyz*b.wxy + a.wxzt*b.wxt + a.wyzt*b.wyt + a.xyzt*b.xyt ; 
c.t =  + a.w*b.wt + a.x*b.xt + a.y*b.yt + a.z*b.zt - a.wt*b.w - a.xt*b.x - a.yt*b.y - a.zt*b.z - a.wxy*b.wxyt - a.wxz*b.wxzt - a.wyz*b.wyzt - a.xyz*b.xyzt + a.wxyt*b.wxy + a.wxzt*b.wxz + a.wyzt*b.wyz + a.xyzt*b.xyz ; 
c.wx =  + a.w*b.x - a.x*b.w - a.wy*b.xy - a.wz*b.xz + a.wt*b.xt + a.xy*b.wy + a.xz*b.wz - a.xt*b.wt - a.wyz*b.xyz + a.wyt*b.xyt + a.wzt*b.xzt + a.xyz*b.wyz - a.xyt*b.wyt - a.xzt*b.wzt - a.wyzt*b.xyzt + a.xyzt*b.wyzt ; 
c.wy =  + a.w*b.y - a.y*b.w + a.wx*b.xy - a.wz*b.yz + a.wt*b.yt - a.xy*b.wx + a.yz*b.wz - a.yt*b.wt + a.wxz*b.xyz - a.wxt*b.xyt + a.wzt*b.yzt - a.xyz*b.wxz + a.xyt*b.wxt - a.yzt*b.wzt + a.wxzt*b.xyzt - a.xyzt*b.wxzt ; 
c.wz =  + a.w*b.z - a.z*b.w + a.wx*b.xz + a.wy*b.yz + a.wt*b.zt - a.xz*b.wx - a.yz*b.wy - a.zt*b.wt - a.wxy*b.xyz - a.wxt*b.xzt - a.wyt*b.yzt + a.xyz*b.wxy + a.xzt*b.wxt + a.yzt*b.wyt - a.wxyt*b.xyzt + a.xyzt*b.wxyt ; 
c.wt =  + a.w*b.t - a.t*b.w + a.wx*b.xt + a.wy*b.yt + a.wz*b.zt - a.xt*b.wx - a.yt*b.wy - a.zt*b.wz - a.wxy*b.xyt - a.wxz*b.xzt - a.wyz*b.yzt + a.xyt*b.wxy + a.xzt*b.wxz + a.yzt*b.wyz - a.wxyz*b.xyzt + a.xyzt*b.wxyz ; 
c.xy =  + a.x*b.y - a.y*b.x - a.wx*b.wy + a.wy*b.wx - a.xz*b.yz + a.xt*b.yt + a.yz*b.xz - a.yt*b.xt - a.wxz*b.wyz + a.wxt*b.wyt + a.wyz*b.wxz - a.wyt*b.wxt + a.xzt*b.yzt - a.yzt*b.xzt - a.wxzt*b.wyzt + a.wyzt*b.wxzt ; 
c.xz =  + a.x*b.z - a.z*b.x - a.wx*b.wz + a.wz*b.wx + a.xy*b.yz + a.xt*b.zt - a.yz*b.xy - a.zt*b.xt + a.wxy*b.wyz + a.wxt*b.wzt - a.wyz*b.wxy - a.wzt*b.wxt - a.xyt*b.yzt + a.yzt*b.xyt + a.wxyt*b.wyzt - a.wyzt*b.wxyt ; 
c.xt =  + a.x*b.t - a.t*b.x - a.wx*b.wt + a.wt*b.wx + a.xy*b.yt + a.xz*b.zt - a.yt*b.xy - a.zt*b.xz + a.wxy*b.wyt + a.wxz*b.wzt - a.wyt*b.wxy - a.wzt*b.wxz - a.xyz*b.yzt + a.yzt*b.xyz + a.wxyz*b.wyzt - a.wyzt*b.wxyz ; 
c.yz =  + a.y*b.z - a.z*b.y - a.wy*b.wz + a.wz*b.wy - a.xy*b.xz + a.xz*b.xy + a.yt*b.zt - a.zt*b.yt - a.wxy*b.wxz + a.wxz*b.wxy + a.wyt*b.wzt - a.wzt*b.wyt + a.xyt*b.xzt - a.xzt*b.xyt - a.wxyt*b.wxzt + a.wxzt*b.wxyt ; 
c.yt =  + a.y*b.t - a.t*b.y - a.wy*b.wt + a.wt*b.wy - a.xy*b.xt + a.xt*b.xy + a.yz*b.zt - a.zt*b.yz - a.wxy*b.wxt + a.wxt*b.wxy + a.wyz*b.wzt - a.wzt*b.wyz + a.xyz*b.xzt - a.xzt*b.xyz - a.wxyz*b.wxzt + a.wxzt*b.wxyz ; 
c.zt =  + a.z*b.t - a.t*b.z - a.wz*b.wt + a.wt*b.wz - a.xz*b.xt + a.xt*b.xz - a.yz*b.yt + a.yt*b.yz - a.wxz*b.wxt + a.wxt*b.wxz - a.wyz*b.wyt + a.wyt*b.wyz - a.xyz*b.xyt + a.xyt*b.xyz + a.wxyz*b.wxyt - a.wxyt*b.wxyz ; 
c.wxy =  - a.z*b.wxyz + a.t*b.wxyt + a.wz*b.xyz - a.wt*b.xyt - a.xz*b.wyz + a.xt*b.wyt + a.yz*b.wxz - a.yt*b.wxt - a.wxz*b.yz + a.wxt*b.yt + a.wyz*b.xz - a.wyt*b.xt - a.xyz*b.wz + a.xyt*b.wt + a.wxyz*b.z - a.wxyt*b.t ; 
c.wxz =  + a.y*b.wxyz + a.t*b.wxzt - a.wy*b.xyz - a.wt*b.xzt + a.xy*b.wyz + a.xt*b.wzt - a.yz*b.wxy - a.zt*b.wxt + a.wxy*b.yz + a.wxt*b.zt - a.wyz*b.xy - a.wzt*b.xt + a.xyz*b.wy + a.xzt*b.wt - a.wxyz*b.y - a.wxzt*b.t ; 
c.wxt =  + a.y*b.wxyt + a.z*b.wxzt - a.wy*b.xyt - a.wz*b.xzt + a.xy*b.wyt + a.xz*b.wzt - a.yt*b.wxy - a.zt*b.wxz + a.wxy*b.yt + a.wxz*b.zt - a.wyt*b.xy - a.wzt*b.xz + a.xyt*b.wy + a.xzt*b.wz - a.wxyt*b.y - a.wxzt*b.z ; 
c.wyz =  - a.x*b.wxyz + a.t*b.wyzt + a.wx*b.xyz - a.wt*b.yzt - a.xy*b.wxz + a.xz*b.wxy + a.yt*b.wzt - a.zt*b.wyt - a.wxy*b.xz + a.wxz*b.xy + a.wyt*b.zt - a.wzt*b.yt - a.xyz*b.wx + a.yzt*b.wt + a.wxyz*b.x - a.wyzt*b.t ; 
c.wyt =  - a.x*b.wxyt + a.z*b.wyzt + a.wx*b.xyt - a.wz*b.yzt - a.xy*b.wxt + a.xt*b.wxy + a.yz*b.wzt - a.zt*b.wyz - a.wxy*b.xt + a.wxt*b.xy + a.wyz*b.zt - a.wzt*b.yz - a.xyt*b.wx + a.yzt*b.wz + a.wxyt*b.x - a.wyzt*b.z ; 
c.wzt =  - a.x*b.wxzt - a.y*b.wyzt + a.wx*b.xzt + a.wy*b.yzt - a.xz*b.wxt + a.xt*b.wxz - a.yz*b.wyt + a.yt*b.wyz - a.wxz*b.xt + a.wxt*b.xz - a.wyz*b.yt + a.wyt*b.yz - a.xzt*b.wx - a.yzt*b.wy + a.wxzt*b.x + a.wyzt*b.y ; 
c.xyz =  + a.w*b.wxyz + a.t*b.xyzt - a.wx*b.wyz + a.wy*b.wxz - a.wz*b.wxy - a.xt*b.yzt + a.yt*b.xzt - a.zt*b.xyt + a.wxy*b.wz - a.wxz*b.wy + a.wyz*b.wx + a.xyt*b.zt - a.xzt*b.yt + a.yzt*b.xt - a.wxyz*b.w - a.xyzt*b.t ; 
c.xyt =  + a.w*b.wxyt + a.z*b.xyzt - a.wx*b.wyt + a.wy*b.wxt - a.wt*b.wxy - a.xz*b.yzt + a.yz*b.xzt - a.zt*b.xyz + a.wxy*b.wt - a.wxt*b.wy + a.wyt*b.wx + a.xyz*b.zt - a.xzt*b.yz + a.yzt*b.xz - a.wxyt*b.w - a.xyzt*b.z ; 
c.xzt =  + a.w*b.wxzt - a.y*b.xyzt - a.wx*b.wzt + a.wz*b.wxt - a.wt*b.wxz + a.xy*b.yzt - a.yz*b.xyt + a.yt*b.xyz + a.wxz*b.wt - a.wxt*b.wz + a.wzt*b.wx - a.xyz*b.yt + a.xyt*b.yz - a.yzt*b.xy - a.wxzt*b.w + a.xyzt*b.y ; 
c.yzt =  + a.w*b.wyzt + a.x*b.xyzt - a.wy*b.wzt + a.wz*b.wyt - a.wt*b.wyz - a.xy*b.xzt + a.xz*b.xyt - a.xt*b.xyz + a.wyz*b.wt - a.wyt*b.wz + a.wzt*b.wy + a.xyz*b.xt - a.xyt*b.xz + a.xzt*b.xy - a.wyzt*b.w - a.xyzt*b.x ; 
c.wxyz =  + a.w*b.xyz - a.x*b.wyz + a.y*b.wxz - a.z*b.wxy + a.wt*b.xyzt - a.xt*b.wyzt + a.yt*b.wxzt - a.zt*b.wxyt + a.wxy*b.z - a.wxz*b.y + a.wyz*b.x - a.xyz*b.w + a.wxyt*b.zt - a.wxzt*b.yt + a.wyzt*b.xt - a.xyzt*b.wt ; 
c.wxyt =  + a.w*b.xyt - a.x*b.wyt + a.y*b.wxt - a.t*b.wxy + a.wz*b.xyzt - a.xz*b.wyzt + a.yz*b.wxzt - a.zt*b.wxyz + a.wxy*b.t - a.wxt*b.y + a.wyt*b.x - a.xyt*b.w + a.wxyz*b.zt - a.wxzt*b.yz + a.wyzt*b.xz - a.xyzt*b.wz ; 
c.wxzt =  + a.w*b.xzt - a.x*b.wzt + a.z*b.wxt - a.t*b.wxz - a.wy*b.xyzt + a.xy*b.wyzt - a.yz*b.wxyt + a.yt*b.wxyz + a.wxz*b.t - a.wxt*b.z + a.wzt*b.x - a.xzt*b.w - a.wxyz*b.yt + a.wxyt*b.yz - a.wyzt*b.xy + a.xyzt*b.wy ; 
c.wyzt =  + a.w*b.yzt - a.y*b.wzt + a.z*b.wyt - a.t*b.wyz + a.wx*b.xyzt - a.xy*b.wxzt + a.xz*b.wxyt - a.xt*b.wxyz + a.wyz*b.t - a.wyt*b.z + a.wzt*b.y - a.yzt*b.w + a.wxyz*b.xt - a.wxyt*b.xz + a.wxzt*b.xy - a.xyzt*b.wx ; 
c.xyzt =  + a.x*b.yzt - a.y*b.xzt + a.z*b.xyt - a.t*b.xyz - a.wx*b.wyzt + a.wy*b.wxzt - a.wz*b.wxyt + a.wt*b.wxyz + a.xyz*b.t - a.xyt*b.z + a.xzt*b.y - a.yzt*b.x - a.wxyz*b.wt + a.wxyt*b.wz - a.wxzt*b.wy + a.wyzt*b.wx ; 
c.wxyzt =  0; 

	return c;
}



//////////////////////////////////////////////////////

GA5_4_1 Inner(GA5_4_1 a, GA5_4_1 b) {

	GA5_4_1 c;

c.q = -a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.t*b.wt - a.wt*b.t + a.wx*b.x + a.wxt*b.xt - a.wxy*b.xy + a.wxyt*b.xyt - a.wxyz*b.xyz - a.wxyzt*b.xyzt - a.wxz*b.xz + a.wxzt*b.xzt + a.wy*b.y + a.wyt*b.yt - a.wyz*b.yz + a.wyzt*b.yzt + a.wz*b.z + a.wzt*b.zt - a.x*b.wx + a.xt*b.wxt - a.xy*b.wxy - a.xyt*b.wxyt + a.xyz*b.wxyz - a.xyzt*b.wxyzt - a.xz*b.wxz - a.xzt*b.wxzt - a.y*b.wy + a.yt*b.wyt - a.yz*b.wyz - a.yzt*b.wyzt - a.z*b.wz + a.zt*b.wzt) ; 

c.x = + (a.t*b.xt + a.w*b.wx - a.wt*b.wxt - a.wx*b.w - a.wxt*b.wt + a.wxy*b.wy - a.wxyt*b.wyt + a.wxyz*b.wyz + a.wxyzt*b.wyzt + a.wxz*b.wz - a.wxzt*b.wzt + a.wy*b.wxy + a.wyt*b.wxyt - a.wyz*b.wxyz + a.wyzt*b.wxyzt + a.wz*b.wxz + a.wzt*b.wxzt - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt) ; 

c.y = + (a.t*b.yt + a.w*b.wy - a.wt*b.wyt - a.wx*b.wxy - a.wxt*b.wxyt - a.wxy*b.wx + a.wxyt*b.wxt - a.wxyz*b.wxz - a.wxyzt*b.wxzt + a.wxz*b.wxyz - a.wxzt*b.wxyzt - a.wy*b.w - a.wyt*b.wt + a.wyz*b.wz - a.wyzt*b.wzt + a.wz*b.wyz + a.wzt*b.wyzt + a.x*b.xy - a.xt*b.xyt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.xz*b.xyz + a.xzt*b.xyzt - a.yt*b.t + a.yz*b.z + a.yzt*b.zt - a.z*b.yz + a.zt*b.yzt) ; 

c.z = + (a.t*b.zt + a.w*b.wz - a.wt*b.wzt - a.wx*b.wxz - a.wxt*b.wxzt - a.wxy*b.wxyz + a.wxyt*b.wxyzt + a.wxyz*b.wxy + a.wxyzt*b.wxyt - a.wxz*b.wx + a.wxzt*b.wxt - a.wy*b.wyz - a.wyt*b.wyzt - a.wyz*b.wy + a.wyzt*b.wyt - a.wz*b.w - a.wzt*b.wt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt + a.y*b.yz - a.yt*b.yzt - a.yz*b.y - a.yzt*b.yt - a.zt*b.t) ; 

c.t = + (a.w*b.wt - a.wt*b.w - a.wx*b.wxt - a.wxt*b.wx - a.wxy*b.wxyt + a.wxyt*b.wxy + a.wxyz*b.wxyzt + a.wxyzt*b.wxyz - a.wxz*b.wxzt + a.wxzt*b.wxz - a.wy*b.wyt - a.wyt*b.wy - a.wyz*b.wyzt + a.wyzt*b.wyz - a.wz*b.wzt - a.wzt*b.wz + a.x*b.xt - a.xt*b.x - a.xy*b.xyt - a.xyt*b.xy - a.xyz*b.xyzt + a.xyzt*b.xyz - a.xz*b.xzt - a.xzt*b.xz + a.y*b.yt - a.yt*b.y - a.yz*b.yzt - a.yzt*b.yz + a.z*b.zt - a.zt*b.z) ; 



c.wx = + (-a.t*b.wxt - a.wxt*b.t + a.wxy*b.y + a.wxyt*b.yt - a.wxyz*b.yz + a.wxyzt*b.yzt + a.wxz*b.z + a.wxzt*b.zt + a.y*b.wxy + a.yt*b.wxyt - a.yz*b.wxyz + a.yzt*b.wxyzt + a.z*b.wxz + a.zt*b.wxzt) ; 

c.wy = + (-a.t*b.wyt - a.wxy*b.x - a.wxyt*b.xt + a.wxyz*b.xz - a.wxyzt*b.xzt - a.wyt*b.t + a.wyz*b.z + a.wyzt*b.zt - a.x*b.wxy - a.xt*b.wxyt + a.xz*b.wxyz - a.xzt*b.wxyzt + a.z*b.wyz + a.zt*b.wyzt) ; 

c.wz = + (-a.t*b.wzt - a.wxyz*b.xy + a.wxyzt*b.xyt - a.wxz*b.x - a.wxzt*b.xt - a.wyz*b.y - a.wyzt*b.yt - a.wzt*b.t - a.x*b.wxz - a.xt*b.wxzt - a.xy*b.wxyz + a.xyt*b.wxyzt - a.y*b.wyz - a.yt*b.wyzt) ; 

c.wt = + (-a.wxt*b.x - a.wxyt*b.xy + a.wxyzt*b.xyz - a.wxzt*b.xz - a.wyt*b.y - a.wyzt*b.yz - a.wzt*b.z - a.x*b.wxt - a.xy*b.wxyt + a.xyz*b.wxyzt - a.xz*b.wxzt - a.y*b.wyt - a.yz*b.wyzt - a.z*b.wzt) ; 

c.xy = + (-a.t*b.xyt + a.w*b.wxy + a.wt*b.wxyt + a.wxy*b.w + a.wxyt*b.wt - a.wxyz*b.wz + a.wxyzt*b.wzt - a.wz*b.wxyz + a.wzt*b.wxyzt - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt + a.z*b.xyz + a.zt*b.xyzt) ; 

c.xz = + (-a.t*b.xzt + a.w*b.wxz + a.wt*b.wxzt + a.wxyz*b.wy - a.wxyzt*b.wyt + a.wxz*b.w + a.wxzt*b.wt + a.wy*b.wxyz - a.wyt*b.wxyzt - a.xyz*b.y - a.xyzt*b.yt - a.xzt*b.t - a.y*b.xyz - a.yt*b.xyzt) ; 

c.xt = + (a.w*b.wxt + a.wxt*b.w + a.wxyt*b.wy - a.wxyzt*b.wyz + a.wxzt*b.wz + a.wy*b.wxyt - a.wyz*b.wxyzt + a.wz*b.wxzt - a.xyt*b.y - a.xyzt*b.yz - a.xzt*b.z - a.y*b.xyt - a.yz*b.xyzt - a.z*b.xzt) ; 

c.yz = + (-a.t*b.yzt + a.w*b.wyz + a.wt*b.wyzt - a.wx*b.wxyz + a.wxt*b.wxyzt - a.wxyz*b.wx + a.wxyzt*b.wxt + a.wyz*b.w + a.wyzt*b.wt + a.x*b.xyz + a.xt*b.xyzt + a.xyz*b.x + a.xyzt*b.xt - a.yzt*b.t) ; 

c.yt = + (a.w*b.wyt - a.wx*b.wxyt - a.wxyt*b.wx + a.wxyzt*b.wxz + a.wxz*b.wxyzt + a.wyt*b.w + a.wyzt*b.wz + a.wz*b.wyzt + a.x*b.xyt + a.xyt*b.x + a.xyzt*b.xz + a.xz*b.xyzt - a.yzt*b.z - a.z*b.yzt) ; 

c.zt = + (a.w*b.wzt - a.wx*b.wxzt - a.wxy*b.wxyzt - a.wxyzt*b.wxy - a.wxzt*b.wx - a.wy*b.wyzt - a.wyzt*b.wy + a.wzt*b.w + a.x*b.xzt - a.xy*b.xyzt - a.xyzt*b.xy + a.xzt*b.x + a.y*b.yzt + a.yzt*b.y) ; 


	c.wxy =  + (a.t*b.wxyt - a.wxyt*b.t + a.wxyz*b.z + a.wxyzt*b.zt - a.z*b.wxyz + a.zt*b.wxyzt) ; 
	c.wxz = + (a.t*b.wxzt - a.wxyz*b.y - a.wxyzt*b.yt - a.wxzt*b.t + a.y*b.wxyz - a.yt*b.wxyzt) ; 
	c.wxt = + (-a.wxyt*b.y - a.wxyzt*b.yz - a.wxzt*b.z + a.y*b.wxyt - a.yz*b.wxyzt + a.z*b.wxzt) ; 
	c.wyz = + (a.t*b.wyzt + a.wxyz*b.x + a.wxyzt*b.xt - a.wyzt*b.t - a.x*b.wxyz + a.xt*b.wxyzt) ; 
	c.wyt = + (a.wxyt*b.x + a.wxyzt*b.xz - a.wyzt*b.z - a.x*b.wxyt + a.xz*b.wxyzt + a.z*b.wyzt) ; 
	c.wzt = + (-a.wxyzt*b.xy + a.wxzt*b.x + a.wyzt*b.y - a.x*b.wxzt - a.xy*b.wxyzt - a.y*b.wyzt) ; 
	c.xyz = + (a.t*b.xyzt + a.w*b.wxyz - a.wt*b.wxyzt - a.wxyz*b.w - a.wxyzt*b.wt - a.xyzt*b.t) ; 
	c.xyt = + (a.w*b.wxyt - a.wxyt*b.w - a.wxyzt*b.wz - a.wz*b.wxyzt - a.xyzt*b.z + a.z*b.xyzt) ; 
	c.xzt = + (a.w*b.wxzt + a.wxyzt*b.wy - a.wxzt*b.w + a.wy*b.wxyzt + a.xyzt*b.y - a.y*b.xyzt) ; 
	c.yzt = + (a.w*b.wyzt - a.wx*b.wxyzt - a.wxyzt*b.wx - a.wyzt*b.w + a.x*b.xyzt - a.xyzt*b.x) ; 

	c.wxyz = + (-a.t*b.wxyzt - a.wxyzt*b.t) ; 
	c.wxyt = + (-a.wxyzt*b.z - a.z*b.wxyzt) ; 
	c.wxzt = + (a.wxyzt*b.y + a.y*b.wxyzt) ; 
	c.wyzt = + (-a.wxyzt*b.x - a.x*b.wxyzt) ; 
	c.xyzt = + (a.w*b.wxyzt + a.wxyzt*b.w) ; 
	
	c.wxyzt = 0 ; 

	return c;
}



//////////////////////////////////////////////////////

GA5_4_1 LeftContraction(GA5_4_1 a, GA5_4_1 b) {

	GA5_4_1 c;

c.q = a.q*b.q - a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.q*b.w + a.t*b.wt - a.x*b.wx + a.xt*b.wxt - a.xy*b.wxy - a.xyt*b.wxyt + a.xyz*b.wxyz - a.xyzt*b.wxyzt - a.xz*b.wxz - a.xzt*b.wxzt - a.y*b.wy + a.yt*b.wyt - a.yz*b.wyz - a.yzt*b.wyzt - a.z*b.wz + a.zt*b.wzt) ; 

c.x = + (a.q*b.x + a.t*b.xt + a.w*b.wx - a.wt*b.wxt + a.wy*b.wxy + a.wyt*b.wxyt - a.wyz*b.wxyz + a.wyzt*b.wxyzt + a.wz*b.wxz + a.wzt*b.wxzt - a.y*b.xy + a.yt*b.xyt - a.yz*b.xyz - a.yzt*b.xyzt - a.z*b.xz + a.zt*b.xzt) ; 

c.y = + (a.q*b.y + a.t*b.yt + a.w*b.wy - a.wt*b.wyt - a.wx*b.wxy - a.wxt*b.wxyt + a.wxz*b.wxyz - a.wxzt*b.wxyzt + a.wz*b.wyz + a.wzt*b.wyzt + a.x*b.xy - a.xt*b.xyt + a.xz*b.xyz + a.xzt*b.xyzt - a.z*b.yz + a.zt*b.yzt) ; 

c.z = + (a.q*b.z + a.t*b.zt + a.w*b.wz - a.wt*b.wzt - a.wx*b.wxz - a.wxt*b.wxzt - a.wxy*b.wxyz + a.wxyt*b.wxyzt - a.wy*b.wyz - a.wyt*b.wyzt + a.x*b.xz - a.xt*b.xzt - a.xy*b.xyz - a.xyt*b.xyzt + a.y*b.yz - a.yt*b.yzt) ; 

c.t = + (a.q*b.t + a.w*b.wt - a.wx*b.wxt - a.wxy*b.wxyt + a.wxyz*b.wxyzt - a.wxz*b.wxzt - a.wy*b.wyt - a.wyz*b.wyzt - a.wz*b.wzt + a.x*b.xt - a.xy*b.xyt - a.xyz*b.xyzt - a.xz*b.xzt + a.y*b.yt - a.yz*b.yzt + a.z*b.zt) ; 

	c.wx = + (a.q*b.wx - a.t*b.wxt + a.y*b.wxy + a.yt*b.wxyt - a.yz*b.wxyz + a.yzt*b.wxyzt + a.z*b.wxz + a.zt*b.wxzt) ; 
	c.wy = + (a.q*b.wy - a.t*b.wyt - a.x*b.wxy - a.xt*b.wxyt + a.xz*b.wxyz - a.xzt*b.wxyzt + a.z*b.wyz + a.zt*b.wyzt) ; 
	c.wz = + (a.q*b.wz - a.t*b.wzt - a.x*b.wxz - a.xt*b.wxzt - a.xy*b.wxyz + a.xyt*b.wxyzt - a.y*b.wyz - a.yt*b.wyzt) ; 
	c.wt = + (a.q*b.wt - a.x*b.wxt - a.xy*b.wxyt + a.xyz*b.wxyzt - a.xz*b.wxzt - a.y*b.wyt - a.yz*b.wyzt - a.z*b.wzt) ; 
	c.xy = + (a.q*b.xy - a.t*b.xyt + a.w*b.wxy + a.wt*b.wxyt - a.wz*b.wxyz + a.wzt*b.wxyzt + a.z*b.xyz + a.zt*b.xyzt) ; 
	c.xz = + (a.q*b.xz - a.t*b.xzt + a.w*b.wxz + a.wt*b.wxzt + a.wy*b.wxyz - a.wyt*b.wxyzt - a.y*b.xyz - a.yt*b.xyzt) ; 
	c.xt = + (a.q*b.xt + a.w*b.wxt + a.wy*b.wxyt - a.wyz*b.wxyzt + a.wz*b.wxzt - a.y*b.xyt - a.yz*b.xyzt - a.z*b.xzt) ; 
	c.yz = + (a.q*b.yz - a.t*b.yzt + a.w*b.wyz + a.wt*b.wyzt - a.wx*b.wxyz + a.wxt*b.wxyzt + a.x*b.xyz + a.xt*b.xyzt) ; 
	c.yt = + (a.q*b.yt + a.w*b.wyt - a.wx*b.wxyt + a.wxz*b.wxyzt + a.wz*b.wyzt + a.x*b.xyt + a.xz*b.xyzt - a.z*b.yzt) ; 
	c.zt = + (a.q*b.zt + a.w*b.wzt - a.wx*b.wxzt - a.wxy*b.wxyzt - a.wy*b.wyzt + a.x*b.xzt - a.xy*b.xyzt + a.y*b.yzt) ; 

	c.wxy = + (a.q*b.wxy + a.t*b.wxyt - a.z*b.wxyz + a.zt*b.wxyzt) ; 
	c.wxz = + (a.q*b.wxz + a.t*b.wxzt + a.y*b.wxyz - a.yt*b.wxyzt) ; 
	c.wxt = + (a.q*b.wxt + a.y*b.wxyt - a.yz*b.wxyzt + a.z*b.wxzt) ; 
	c.wyz = + (a.q*b.wyz + a.t*b.wyzt - a.x*b.wxyz + a.xt*b.wxyzt) ; 
	c.wyt = + (a.q*b.wyt - a.x*b.wxyt + a.xz*b.wxyzt + a.z*b.wyzt) ; 
	c.wzt = + (a.q*b.wzt - a.x*b.wxzt - a.xy*b.wxyzt - a.y*b.wyzt) ; 
	c.xyz = + (a.q*b.xyz + a.t*b.xyzt + a.w*b.wxyz - a.wt*b.wxyzt) ; 
	c.xyt = + (a.q*b.xyt + a.w*b.wxyt - a.wz*b.wxyzt + a.z*b.xyzt) ; 
	c.xzt = + (a.q*b.xzt + a.w*b.wxzt + a.wy*b.wxyzt - a.y*b.xyzt) ; 
	c.yzt = + (a.q*b.yzt + a.w*b.wyzt - a.wx*b.wxyzt + a.x*b.xyzt) ; 

	c.wxyz = + (a.q*b.wxyz - a.t*b.wxyzt) ; 
	c.wxyt = + (a.q*b.wxyt - a.z*b.wxyzt) ; 
	c.wxzt = + (a.q*b.wxzt + a.y*b.wxyzt) ; 
	c.wyzt = + (a.q*b.wyzt - a.x*b.wxyzt) ; 
	c.xyzt = + (a.q*b.xyzt + a.w*b.wxyzt) ; 

	c.wxyzt = + a.q*b.wxyzt ; 

	return c;
}


//////////////////////////////////////////////////////

GA5_4_1 RightContraction(GA5_4_1 a, GA5_4_1 b) {

	GA5_4_1 c;

c.q = a.q*b.q - a.t*b.t + a.w*b.w + a.wt*b.wt - a.wx*b.wx + a.wxt*b.wxt - a.wxy*b.wxy - a.wxyt*b.wxyt + a.wxyz*b.wxyz - a.wxyzt*b.wxyzt - a.wxz*b.wxz - a.wxzt*b.wxzt - a.wy*b.wy + a.wyt*b.wyt - a.wyz*b.wyz - a.wyzt*b.wyzt - a.wz*b.wz + a.wzt*b.wzt + a.x*b.x + a.xt*b.xt - a.xy*b.xy + a.xyt*b.xyt - a.xyz*b.xyz - a.xyzt*b.xyzt - a.xz*b.xz + a.xzt*b.xzt + a.y*b.y + a.yt*b.yt - a.yz*b.yz + a.yzt*b.yzt + a.z*b.z + a.zt*b.zt ; 

c.w = + (a.w*b.q - a.wt*b.t + a.wx*b.x + a.wxt*b.xt - a.wxy*b.xy + a.wxyt*b.xyt - a.wxyz*b.xyz - a.wxyzt*b.xyzt - a.wxz*b.xz + a.wxzt*b.xzt + a.wy*b.y + a.wyt*b.yt - a.wyz*b.yz + a.wyzt*b.yzt + a.wz*b.z + a.wzt*b.zt) ; 

c.x = + (-a.wx*b.w - a.wxt*b.wt + a.wxy*b.wy - a.wxyt*b.wyt + a.wxyz*b.wyz + a.wxyzt*b.wyzt + a.wxz*b.wz - a.wxzt*b.wzt + a.x*b.q - a.xt*b.t + a.xy*b.y + a.xyt*b.yt - a.xyz*b.yz + a.xyzt*b.yzt + a.xz*b.z + a.xzt*b.zt) ; 

c.y = + (-a.wxy*b.wx + a.wxyt*b.wxt - a.wxyz*b.wxz - a.wxyzt*b.wxzt - a.wy*b.w - a.wyt*b.wt + a.wyz*b.wz - a.wyzt*b.wzt - a.xy*b.x - a.xyt*b.xt + a.xyz*b.xz - a.xyzt*b.xzt + a.y*b.q - a.yt*b.t + a.yz*b.z + a.yzt*b.zt) ; 

c.z = + (a.wxyz*b.wxy + a.wxyzt*b.wxyt - a.wxz*b.wx + a.wxzt*b.wxt - a.wyz*b.wy + a.wyzt*b.wyt - a.wz*b.w - a.wzt*b.wt - a.xyz*b.xy + a.xyzt*b.xyt - a.xz*b.x - a.xzt*b.xt - a.yz*b.y - a.yzt*b.yt + a.z*b.q - a.zt*b.t) ; 

c.t = + (a.t*b.q - a.wt*b.w - a.wxt*b.wx + a.wxyt*b.wxy + a.wxyzt*b.wxyz + a.wxzt*b.wxz - a.wyt*b.wy + a.wyzt*b.wyz - a.wzt*b.wz - a.xt*b.x - a.xyt*b.xy + a.xyzt*b.xyz - a.xzt*b.xz - a.yt*b.y - a.yzt*b.yz - a.zt*b.z) ; 

	c.wx = + (a.wx*b.q - a.wxt*b.t + a.wxy*b.y + a.wxyt*b.yt - a.wxyz*b.yz + a.wxyzt*b.yzt + a.wxz*b.z + a.wxzt*b.zt) ; 
	c.wy = + (-a.wxy*b.x - a.wxyt*b.xt + a.wxyz*b.xz - a.wxyzt*b.xzt + a.wy*b.q - a.wyt*b.t + a.wyz*b.z + a.wyzt*b.zt) ; 
	c.wz = + (-a.wxyz*b.xy + a.wxyzt*b.xyt - a.wxz*b.x - a.wxzt*b.xt - a.wyz*b.y - a.wyzt*b.yt + a.wz*b.q - a.wzt*b.t) ; 
	c.wt = + (a.wt*b.q - a.wxt*b.x - a.wxyt*b.xy + a.wxyzt*b.xyz - a.wxzt*b.xz - a.wyt*b.y - a.wyzt*b.yz - a.wzt*b.z) ; 
	c.xy = + (a.wxy*b.w + a.wxyt*b.wt - a.wxyz*b.wz + a.wxyzt*b.wzt + a.xy*b.q - a.xyt*b.t + a.xyz*b.z + a.xyzt*b.zt) ; 
	c.xz = + (a.wxyz*b.wy - a.wxyzt*b.wyt + a.wxz*b.w + a.wxzt*b.wt - a.xyz*b.y - a.xyzt*b.yt + a.xz*b.q - a.xzt*b.t) ; 
	c.xt = + (a.wxt*b.w + a.wxyt*b.wy - a.wxyzt*b.wyz + a.wxzt*b.wz + a.xt*b.q - a.xyt*b.y - a.xyzt*b.yz - a.xzt*b.z) ; 
	c.yz = + (-a.wxyz*b.wx + a.wxyzt*b.wxt + a.wyz*b.w + a.wyzt*b.wt + a.xyz*b.x + a.xyzt*b.xt + a.yz*b.q - a.yzt*b.t) ; 
	c.yt = + (-a.wxyt*b.wx + a.wxyzt*b.wxz + a.wyt*b.w + a.wyzt*b.wz + a.xyt*b.x + a.xyzt*b.xz + a.yt*b.q - a.yzt*b.z) ; 
	c.zt = + (-a.wxyzt*b.wxy - a.wxzt*b.wx - a.wyzt*b.wy + a.wzt*b.w - a.xyzt*b.xy + a.xzt*b.x + a.yzt*b.y + a.zt*b.q) ; 

	c.wxy = + (a.wxy*b.q - a.wxyt*b.t + a.wxyz*b.z + a.wxyzt*b.zt) ; 
	c.wxz = + (-a.wxyz*b.y - a.wxyzt*b.yt + a.wxz*b.q - a.wxzt*b.t) ; 
	c.wxt = + (a.wxt*b.q - a.wxyt*b.y - a.wxyzt*b.yz - a.wxzt*b.z) ; 
	c.wyz = + (a.wxyz*b.x + a.wxyzt*b.xt + a.wyz*b.q - a.wyzt*b.t) ; 
	c.wyt = + (a.wxyt*b.x + a.wxyzt*b.xz + a.wyt*b.q - a.wyzt*b.z) ; 
	c.wzt = + (-a.wxyzt*b.xy + a.wxzt*b.x + a.wyzt*b.y + a.wzt*b.q) ; 
	c.xyz = + (-a.wxyz*b.w - a.wxyzt*b.wt + a.xyz*b.q - a.xyzt*b.t) ; 
	c.xyt = + (-a.wxyt*b.w - a.wxyzt*b.wz + a.xyt*b.q - a.xyzt*b.z) ; 
	c.xzt = + (a.wxyzt*b.wy - a.wxzt*b.w + a.xyzt*b.y + a.xzt*b.q) ; 
	c.yzt = + (-a.wxyzt*b.wx - a.wyzt*b.w - a.xyzt*b.x + a.yzt*b.q) ; 

	c.wxyz = + (a.wxyz*b.q - a.wxyzt*b.t) ; 
	c.wxyt = + (a.wxyt*b.q - a.wxyzt*b.z) ; 
	c.wxzt = + (a.wxyzt*b.y + a.wxzt*b.q) ; 
	c.wyzt = + (-a.wxyzt*b.x + a.wyzt*b.q) ; 
	c.xyzt = + (a.wxyzt*b.w + a.xyzt*b.q) ; 

	c.wxyzt = + a.wxyzt*b.q ;

	return c;
}

//////////////////////////////////////////////////////

GA5_4_1 Parity(GA5_4_1 a)	// Corresponds to parity transform: w -> -w, x -> -x, y -> -y, z ->-z, t -> -t
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = a.wx;		b.wy = a.wy;		b.wz = a.wz;		b.wt = a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = - a.wxy;	b.wxz = - a.wxz;	b.wxt = - a.wxt;	b.wyz = - a.wyz;	b.wyt = - a.wyt;
	b.wzt = - a.wzt;	b.xyz = - a.xyz;	b.xyt = - a.xyt;	b.xzt = - a.xzt;	b.yzt = - a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt = - a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 ComplexConjugate(GA5_4_1 a)	// negate every term including w
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = - a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = a.t;
	b.wx = - a.wx;		b.wy = - a.wy;		b.wz = - a.wz;		b.wt = - a.wt;		b.xy = a.xy;
	b.xz = a.xz;		b.yz = a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = - a.wxy;	b.wxz = - a.wxz;	b.wxt = - a.wxt;	b.wyz = - a.wyz;	b.wyt = - a.wyt;
	b.wzt = - a.wzt;	b.xyz = a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz = - a.wxyz;	b.wxyt = - a.wxyt;	b.wxzt = - a.wxzt;	b.wyzt = - a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt = - a.wxyzt;

	return b;
}
//////////////////////////////////////////////////////

GA5_4_1 Hermitian(GA5_4_1 a)	
//	+ a*q
//	+ b*w + c*x + d*y + e*z - f*t
//	- g*wx - h*wy - j*wz + k*wt - l*xy - m*xz + n*xt - p*yz + r*yt + s*zt
//	- S*wxy - R*wxz + P*wxt - N*wyz + M*wyt + L*wzt - K*xyz + J*xyt + H*xzt + G*yzt
//	+ F*wxyz - E*wxyt - D*wxzt - C*wyzt - B*xyzt
//	- A*wxyzt 
{
	GA5_4_1 b;

	b.q =  a.q;
	b.w = a.w;		b.x = a.x;		b.y = a.y;		b.z = a.z;		b.t = -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = a.wxt;		b.wyz = -a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = -a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = -a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 Conjugate(GA5_4_1 a)	
{
	GA5_4_1 b;

/*	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = -a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = -a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = -a.xt;		b.yt = -a.yt;		b.zt = -a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = -a.wxt;		b.wyz = -a.wyz;		b.wyt = -a.wyt;
	b.wzt = -a.wzt;		b.xyz = -a.xyz;		b.xyt = -a.xyt;		b.xzt = -a.xzt;		b.yzt = -a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  -a.wxyt;	b.wxzt =  -a.wxzt;	b.wyzt =  -a.wyzt;	b.xyzt =  -a.xyzt;
	b.wxyzt = -a.wxyzt;
*/
	b.q =  a.q;
	b.w = -a.w;		b.x = -a.x;		b.y = -a.y;		b.z = -a.z;		b.t = a.t;
	b.wx = -a.wx;		b.wy = -a.wy;		b.wz = -a.wz;		b.wt = a.wt;		b.xy = -a.xy;
	b.xz = -a.xz;		b.yz = -a.yz;		b.xt = a.xt;		b.yt = a.yt;		b.zt = a.zt;
	b.wxy = -a.wxy;		b.wxz = -a.wxz;		b.wxt = a.wxt;		b.wyz = -a.wyz;		b.wyt = a.wyt;
	b.wzt = a.wzt;		b.xyz = -a.xyz;		b.xyt = a.xyt;		b.xzt = a.xzt;		b.yzt = a.yzt;
	b.wxyz =  -a.wxyz;	b.wxyt =  a.wxyt;	b.wxzt =  a.wxzt;	b.wyzt =  a.wyzt;	b.xyzt =  a.xyzt;
	b.wxyzt = a.wxyzt;

	return b;
}

//////////////////////////////////////////////////////

GA5_4_1 SymmetricProductViaFormula(GA5_4_1 u, GA5_4_1 v)	// symmetric product
{
	GA5_4_1 a;
	a = Divide_By_Constant(Add(Product(u,v) , Product(v,u)),2.0);
	return a;
}

//////////////////////////////////////////////////////

GA5_4_1 AntiSymmetricProductViaFormula(GA5_4_1 u, GA5_4_1 v)	// antisymmetric product
{
	GA5_4_1 a;
	a = Divide_By_Constant(Subtract(Product(u,v) , Product(v,u)),2.0);
	return a;
}

//////////////////////////////////////////////////////

complex Determinant(GA5_4_1 MV1)
{
	complex det;
	long a, b,c,d,e,f, g,G, F,E,D,C,B, A;
//	GA5_4_1 MV2, MV3, MV4, MV5;

//	MV2 = Reverse(MV1);
//	MV3 = MV1*MV2;	// result has zero in bivector and trivector components

//	MV4 = Reverse_Vector_Quad(MV3);
//	MV5 = MV3*MV4;
//	det = MV5.q + I*MV5.wxyzt;

//	deltat for 100X MV determinant = 10 => 100 msec/determinant using above
//			3.8 msec using code below..   25X improvement


	a = + MV1.q*MV1.q + MV1.w*MV1.w + MV1.x*MV1.x + MV1.y*MV1.y + MV1.z*MV1.z - MV1.t*MV1.t
		 + MV1.wx*MV1.wx + MV1.wy*MV1.wy + MV1.wz*MV1.wz - MV1.wt*MV1.wt + MV1.xy*MV1.xy
		 + MV1.xz*MV1.xz - MV1.xt*MV1.xt + MV1.yz*MV1.yz - MV1.yt*MV1.yt - MV1.zt*MV1.zt
		 + MV1.wxy*MV1.wxy + MV1.wxz*MV1.wxz - MV1.wxt*MV1.wxt + MV1.wyz*MV1.wyz - MV1.wyt*MV1.wyt 
		 - MV1.wzt*MV1.wzt + MV1.xyz*MV1.xyz - MV1.xyt*MV1.xyt - MV1.xzt*MV1.xzt - MV1.yzt*MV1.yzt
		 + MV1.wxyz*MV1.wxyz - MV1.wxyt*MV1.wxyt - MV1.wxzt*MV1.wxzt - MV1.wyzt*MV1.wyzt - MV1.xyzt*MV1.xyzt
		 - MV1.wxyzt*MV1.wxyzt ; 

	b = 2*MV1.y*MV1.wy-2*MV1.xzt*MV1.wxzt+2*MV1.wyz*MV1.yz+2*MV1.wxz*MV1.xz-2*MV1.yzt*MV1.wyzt
		-2*MV1.t*MV1.wt+2*MV1.wxy*MV1.xy-2*MV1.wxyt*MV1.xyt-2*MV1.yt*MV1.wyt-2*MV1.wxt*MV1.xt+2*MV1.z*MV1.wz
		-2*MV1.xyzt*MV1.wxyzt+2*MV1.xyz*MV1.wxyz+2*MV1.x*MV1.wx+2*MV1.w*MV1.q-2*MV1.wzt*MV1.zt ; 

	c = -2*MV1.wyz*MV1.wxyz-2*MV1.xyzt*MV1.yzt-2*MV1.xzt*MV1.zt-2*MV1.w*MV1.wx+2*MV1.wxyzt*MV1.wyzt
		-2*MV1.wxy*MV1.wy+2*MV1.x*MV1.q+2*MV1.wxyt*MV1.wyt+2*MV1.y*MV1.xy+2*MV1.z*MV1.xz+2*MV1.wxt*MV1.wt
		-2*MV1.yt*MV1.xyt-2*MV1.t*MV1.xt-2*MV1.wxz*MV1.wz+2*MV1.xyz*MV1.yz+2*MV1.wzt*MV1.wxzt ; 

	d = 2*MV1.xt*MV1.xyt-2*MV1.wxyt*MV1.wxt-2*MV1.x*MV1.xy+2*MV1.wzt*MV1.wyzt+2*MV1.xyzt*MV1.xzt
		-2*MV1.wz*MV1.wyz+2*MV1.y*MV1.q-2*MV1.w*MV1.wy+2*MV1.wt*MV1.wyt-2*MV1.wxyzt*MV1.wxzt-2*MV1.zt*MV1.yzt
		-2*MV1.xz*MV1.xyz+2*MV1.z*MV1.yz+2*MV1.wxz*MV1.wxyz+2*MV1.wx*MV1.wxy-2*MV1.yt*MV1.t ; 

	e = 2*MV1.wxyt*MV1.wxyzt+2*MV1.xzt*MV1.xt+2*MV1.wzt*MV1.wt-2*MV1.wyzt*MV1.wyt
		+2*MV1.wyz*MV1.wy-2*MV1.xz*MV1.x+2*MV1.yt*MV1.yzt+2*MV1.wxz*MV1.wx-2*MV1.y*MV1.yz-2*MV1.wxy*MV1.wxyz
		-2*MV1.wxt*MV1.wxzt-2*MV1.t*MV1.zt+2*MV1.z*MV1.q-2*MV1.xyzt*MV1.xyt+2*MV1.xyz*MV1.xy-2*MV1.w*MV1.wz ; 

	f = 2*MV1.wy*MV1.wyt+2*MV1.wxt*MV1.wx-2*MV1.xyzt*MV1.xyz-2*MV1.yt*MV1.y-2*MV1.w*MV1.wt
		-2*MV1.x*MV1.xt+2*MV1.xy*MV1.xyt+2*MV1.yzt*MV1.yz+2*MV1.wz*MV1.wzt-2*MV1.wxz*MV1.wxzt-2*MV1.wxyt*MV1.wxy
		-2*MV1.wyz*MV1.wyzt+2*MV1.xz*MV1.xzt+2*MV1.wxyzt*MV1.wxyz-2*MV1.z*MV1.zt+2*MV1.t*MV1.q ; 

	F = -2*MV1.xzt*MV1.wyt+2*MV1.xyzt*MV1.wt-2*MV1.wx*MV1.yz-2*MV1.w*MV1.xyz+2*MV1.wxt*MV1.yzt
		+2*MV1.z*MV1.wxy+2*MV1.xz*MV1.wy-2*MV1.wxyt*MV1.zt-2*MV1.wxz*MV1.y-2*MV1.wz*MV1.xy+2*MV1.yt*MV1.wxzt
		-2*MV1.xt*MV1.wyzt+2*MV1.q*MV1.wxyz-2*MV1.t*MV1.wxyzt+2*MV1.wzt*MV1.xyt+2*MV1.x*MV1.wyz ; 

	E = 2*MV1.wzt*MV1.xyz-2*MV1.yt*MV1.wx+2*MV1.wxzt*MV1.yz+2*MV1.t*MV1.wxy+2*MV1.wxz*MV1.yzt
		+2*MV1.xt*MV1.wy+2*MV1.wxyt*MV1.q-2*MV1.zt*MV1.wxyz-2*MV1.wxt*MV1.y+2*MV1.x*MV1.wyt+2*MV1.xyzt*MV1.wz
		-2*MV1.z*MV1.wxyzt-2*MV1.w*MV1.xyt-2*MV1.xzt*MV1.wyz-2*MV1.wt*MV1.xy-2*MV1.xz*MV1.wyzt ; 

	D = -2*MV1.w*MV1.xzt+2*MV1.wxzt*MV1.q-2*MV1.wxyt*MV1.yz+2*MV1.wxyzt*MV1.y-2*MV1.yzt*MV1.wxy
		+2*MV1.wyz*MV1.xyt+2*MV1.wyzt*MV1.xy-2*MV1.xyz*MV1.wyt-2*MV1.z*MV1.wxt-2*MV1.xz*MV1.wt+2*MV1.wxz*MV1.t
		+2*MV1.wz*MV1.xt-2*MV1.zt*MV1.wx+2*MV1.wzt*MV1.x-2*MV1.xyzt*MV1.wy+2*MV1.yt*MV1.wxyz ; 

	C = 2*MV1.t*MV1.wyz-2*MV1.wxzt*MV1.xy+2*MV1.wzt*MV1.y+2*MV1.xyzt*MV1.wx+2*MV1.wyzt*MV1.q
		-2*MV1.w*MV1.yzt-2*MV1.xt*MV1.wxyz+2*MV1.wxt*MV1.xyz-2*MV1.zt*MV1.wy-2*MV1.wxyzt*MV1.x+2*MV1.xzt*MV1.wxy
		-2*MV1.z*MV1.wyt+2*MV1.yt*MV1.wz-2*MV1.wxz*MV1.xyt+2*MV1.wxyt*MV1.xz-2*MV1.wt*MV1.yz ; 

	B = 2*MV1.yt*MV1.xz-2*MV1.wxt*MV1.wyz-2*MV1.zt*MV1.xy-2*MV1.wxyt*MV1.wz+2*MV1.t*MV1.xyz
		-2*MV1.xt*MV1.yz-2*MV1.wzt*MV1.wxy+2*MV1.wxzt*MV1.wy+2*MV1.xzt*MV1.y+2*MV1.xyzt*MV1.q-2*MV1.x*MV1.yzt
		-2*MV1.wx*MV1.wyzt+2*MV1.w*MV1.wxyzt+2*MV1.wt*MV1.wxyz-2*MV1.z*MV1.xyt+2*MV1.wxz*MV1.wyt ; 

	A = 2*MV1.xzt*MV1.wy+2*MV1.wxz*MV1.yt+2*MV1.xyzt*MV1.w+2*MV1.wt*MV1.xyz+2*MV1.xz*MV1.wyt
		-2*MV1.z*MV1.wxyt-2*MV1.yzt*MV1.wx-2*MV1.xt*MV1.wyz-2*MV1.wzt*MV1.xy-2*MV1.wxt*MV1.yz+2*MV1.t*MV1.wxyz
		-2*MV1.x*MV1.wyzt+2*MV1.wxyzt*MV1.q-2*MV1.zt*MV1.wxy+2*MV1.wxzt*MV1.y-2*MV1.wz*MV1.xyt ; 


//	g = C^2-c^2-F^2-A^2+f^2+a^2+D^2-d^2+B^2-b^2+E^2-e^2 ; 

	g = a*a - b*b -c*c - d*d - e*e + f*f - F*F + E*E + D*D + C*C + B*B - A*A;

	G = + 2*a*A - 2*b*B + 2*c*C - 2*d*D + 2*e*E - 2*f*F ; 

//	det = complex(g, G);
	det.r = g;
	det.i = G;
	
	return det;
}

/////////////////////////////////////////////////////////////////////

GA5_4_1 Adjugate(GA5_4_1 r)
{

// Adjugate = Reverse(r)*Reverse_Vector_Quad(r*Reverse(r));

	long a, b,c,d,e,f, F,E,D,C,B, A;
	GA5_4_1 MV1, MV2, Adj;
	
	MV1 = Reverse(r);
//	MV2 = Reverse_Vector_Quad(r*Reverse(r););
//	Adj = MV1*MV2;
	
	a = + r.q*r.q + r.w*r.w + r.x*r.x + r.y*r.y + r.z*r.z - r.t*r.t
		 + r.wx*r.wx + r.wy*r.wy + r.wz*r.wz - r.wt*r.wt + r.xy*r.xy
		 + r.xz*r.xz - r.xt*r.xt + r.yz*r.yz - r.yt*r.yt - r.zt*r.zt
		 + r.wxy*r.wxy + r.wxz*r.wxz - r.wxt*r.wxt + r.wyz*r.wyz - r.wyt*r.wyt 
		 - r.wzt*r.wzt + r.xyz*r.xyz - r.xyt*r.xyt - r.xzt*r.xzt - r.yzt*r.yzt
		 + r.wxyz*r.wxyz - r.wxyt*r.wxyt - r.wxzt*r.wxzt - r.wyzt*r.wyzt - r.xyzt*r.xyzt
		 - r.wxyzt*r.wxyzt ; 

	b = 2*r.y*r.wy-2*r.xzt*r.wxzt+2*r.wyz*r.yz+2*r.wxz*r.xz-2*r.yzt*r.wyzt
		-2*r.t*r.wt+2*r.wxy*r.xy-2*r.wxyt*r.xyt-2*r.yt*r.wyt-2*r.wxt*r.xt+2*r.z*r.wz
		-2*r.xyzt*r.wxyzt+2*r.xyz*r.wxyz+2*r.x*r.wx+2*r.w*r.q-2*r.wzt*r.zt ; 

	c = -2*r.wyz*r.wxyz-2*r.xyzt*r.yzt-2*r.xzt*r.zt-2*r.w*r.wx+2*r.wxyzt*r.wyzt
		-2*r.wxy*r.wy+2*r.x*r.q+2*r.wxyt*r.wyt+2*r.y*r.xy+2*r.z*r.xz+2*r.wxt*r.wt
		-2*r.yt*r.xyt-2*r.t*r.xt-2*r.wxz*r.wz+2*r.xyz*r.yz+2*r.wzt*r.wxzt ; 

	d = 2*r.xt*r.xyt-2*r.wxyt*r.wxt-2*r.x*r.xy+2*r.wzt*r.wyzt+2*r.xyzt*r.xzt
		-2*r.wz*r.wyz+2*r.y*r.q-2*r.w*r.wy+2*r.wt*r.wyt-2*r.wxyzt*r.wxzt-2*r.zt*r.yzt
		-2*r.xz*r.xyz+2*r.z*r.yz+2*r.wxz*r.wxyz+2*r.wx*r.wxy-2*r.yt*r.t ; 

	e = 2*r.wxyt*r.wxyzt+2*r.xzt*r.xt+2*r.wzt*r.wt-2*r.wyzt*r.wyt
		+2*r.wyz*r.wy-2*r.xz*r.x+2*r.yt*r.yzt+2*r.wxz*r.wx-2*r.y*r.yz-2*r.wxy*r.wxyz
		-2*r.wxt*r.wxzt-2*r.t*r.zt+2*r.z*r.q-2*r.xyzt*r.xyt+2*r.xyz*r.xy-2*r.w*r.wz ; 

	f = 2*r.wy*r.wyt+2*r.wxt*r.wx-2*r.xyzt*r.xyz-2*r.yt*r.y-2*r.w*r.wt
		-2*r.x*r.xt+2*r.xy*r.xyt+2*r.yzt*r.yz+2*r.wz*r.wzt-2*r.wxz*r.wxzt-2*r.wxyt*r.wxy
		-2*r.wyz*r.wyzt+2*r.xz*r.xzt+2*r.wxyzt*r.wxyz-2*r.z*r.zt+2*r.t*r.q ; 

	F = -2*r.xzt*r.wyt+2*r.xyzt*r.wt-2*r.wx*r.yz-2*r.w*r.xyz+2*r.wxt*r.yzt
		+2*r.z*r.wxy+2*r.xz*r.wy-2*r.wxyt*r.zt-2*r.wxz*r.y-2*r.wz*r.xy+2*r.yt*r.wxzt
		-2*r.xt*r.wyzt+2*r.q*r.wxyz-2*r.t*r.wxyzt+2*r.wzt*r.xyt+2*r.x*r.wyz ; 

	E = 2*r.wzt*r.xyz-2*r.yt*r.wx+2*r.wxzt*r.yz+2*r.t*r.wxy+2*r.wxz*r.yzt
		+2*r.xt*r.wy+2*r.wxyt*r.q-2*r.zt*r.wxyz-2*r.wxt*r.y+2*r.x*r.wyt+2*r.xyzt*r.wz
		-2*r.z*r.wxyzt-2*r.w*r.xyt-2*r.xzt*r.wyz-2*r.wt*r.xy-2*r.xz*r.wyzt ; 

	D = -2*r.w*r.xzt+2*r.wxzt*r.q-2*r.wxyt*r.yz+2*r.wxyzt*r.y-2*r.yzt*r.wxy
		+2*r.wyz*r.xyt+2*r.wyzt*r.xy-2*r.xyz*r.wyt-2*r.z*r.wxt-2*r.xz*r.wt+2*r.wxz*r.t
		+2*r.wz*r.xt-2*r.zt*r.wx+2*r.wzt*r.x-2*r.xyzt*r.wy+2*r.yt*r.wxyz ; 

	C = 2*r.t*r.wyz-2*r.wxzt*r.xy+2*r.wzt*r.y+2*r.xyzt*r.wx+2*r.wyzt*r.q
		-2*r.w*r.yzt-2*r.xt*r.wxyz+2*r.wxt*r.xyz-2*r.zt*r.wy-2*r.wxyzt*r.x+2*r.xzt*r.wxy
		-2*r.z*r.wyt+2*r.yt*r.wz-2*r.wxz*r.xyt+2*r.wxyt*r.xz-2*r.wt*r.yz ; 

	B = 2*r.yt*r.xz-2*r.wxt*r.wyz-2*r.zt*r.xy-2*r.wxyt*r.wz+2*r.t*r.xyz
		-2*r.xt*r.yz-2*r.wzt*r.wxy+2*r.wxzt*r.wy+2*r.xzt*r.y+2*r.xyzt*r.q-2*r.x*r.yzt
		-2*r.wx*r.wyzt+2*r.w*r.wxyzt+2*r.wt*r.wxyz-2*r.z*r.xyt+2*r.wxz*r.wyt ; 

	A = 2*r.xzt*r.wy+2*r.wxz*r.yt+2*r.xyzt*r.w+2*r.wt*r.xyz+2*r.xz*r.wyt
		-2*r.z*r.wxyt-2*r.yzt*r.wx-2*r.xt*r.wyz-2*r.wzt*r.xy-2*r.wxt*r.yz+2*r.t*r.wxyz
		-2*r.x*r.wyzt+2*r.wxyzt*r.q-2*r.zt*r.wxy+2*r.wxzt*r.y-2*r.wz*r.xyt ; 

	// MV2 = Reverse_Vector_Quad(r*Reverse(r));
	MV2 = Set(a, -b,-c,-d,-e,-f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, -F,-E,-D,-C,-B, A);
	
	Adj = Product(MV1,MV2);  // Adjugate = Reverse(r)*Reverse_Vector_Quad(r*Reverse(r));	

	return Adj;
}

/////////////////////////////////////////////////////////////////////

GA5_4_1 Reciprocal( GA5_4_1 r)
{
	long det2;		// det squared
	long a, b,c,d,e,f, g,G, F,E,D,C,B, A;
	GA5_4_1 MV1, MV2, MV3, Result;
	
	MV1 = Reverse(r);
//	MV2 = r*Reverse(r);
//	MV3 = Reverse_Vector_Quad(MV2);
	
	a = + r.q*r.q + r.w*r.w + r.x*r.x + r.y*r.y + r.z*r.z - r.t*r.t
		 + r.wx*r.wx + r.wy*r.wy + r.wz*r.wz - r.wt*r.wt + r.xy*r.xy
		 + r.xz*r.xz - r.xt*r.xt + r.yz*r.yz - r.yt*r.yt - r.zt*r.zt
		 + r.wxy*r.wxy + r.wxz*r.wxz - r.wxt*r.wxt + r.wyz*r.wyz - r.wyt*r.wyt 
		 - r.wzt*r.wzt + r.xyz*r.xyz - r.xyt*r.xyt - r.xzt*r.xzt - r.yzt*r.yzt
		 + r.wxyz*r.wxyz - r.wxyt*r.wxyt - r.wxzt*r.wxzt - r.wyzt*r.wyzt - r.xyzt*r.xyzt
		 - r.wxyzt*r.wxyzt ; 

	b = 2*r.y*r.wy-2*r.xzt*r.wxzt+2*r.wyz*r.yz+2*r.wxz*r.xz-2*r.yzt*r.wyzt
		-2*r.t*r.wt+2*r.wxy*r.xy-2*r.wxyt*r.xyt-2*r.yt*r.wyt-2*r.wxt*r.xt+2*r.z*r.wz
		-2*r.xyzt*r.wxyzt+2*r.xyz*r.wxyz+2*r.x*r.wx+2*r.w*r.q-2*r.wzt*r.zt ; 

	c = -2*r.wyz*r.wxyz-2*r.xyzt*r.yzt-2*r.xzt*r.zt-2*r.w*r.wx+2*r.wxyzt*r.wyzt
		-2*r.wxy*r.wy+2*r.x*r.q+2*r.wxyt*r.wyt+2*r.y*r.xy+2*r.z*r.xz+2*r.wxt*r.wt
		-2*r.yt*r.xyt-2*r.t*r.xt-2*r.wxz*r.wz+2*r.xyz*r.yz+2*r.wzt*r.wxzt ; 

	d = 2*r.xt*r.xyt-2*r.wxyt*r.wxt-2*r.x*r.xy+2*r.wzt*r.wyzt+2*r.xyzt*r.xzt
		-2*r.wz*r.wyz+2*r.y*r.q-2*r.w*r.wy+2*r.wt*r.wyt-2*r.wxyzt*r.wxzt-2*r.zt*r.yzt
		-2*r.xz*r.xyz+2*r.z*r.yz+2*r.wxz*r.wxyz+2*r.wx*r.wxy-2*r.yt*r.t ; 

	e = 2*r.wxyt*r.wxyzt+2*r.xzt*r.xt+2*r.wzt*r.wt-2*r.wyzt*r.wyt
		+2*r.wyz*r.wy-2*r.xz*r.x+2*r.yt*r.yzt+2*r.wxz*r.wx-2*r.y*r.yz-2*r.wxy*r.wxyz
		-2*r.wxt*r.wxzt-2*r.t*r.zt+2*r.z*r.q-2*r.xyzt*r.xyt+2*r.xyz*r.xy-2*r.w*r.wz ; 

	f = 2*r.wy*r.wyt+2*r.wxt*r.wx-2*r.xyzt*r.xyz-2*r.yt*r.y-2*r.w*r.wt
		-2*r.x*r.xt+2*r.xy*r.xyt+2*r.yzt*r.yz+2*r.wz*r.wzt-2*r.wxz*r.wxzt-2*r.wxyt*r.wxy
		-2*r.wyz*r.wyzt+2*r.xz*r.xzt+2*r.wxyzt*r.wxyz-2*r.z*r.zt+2*r.t*r.q ; 

	F = -2*r.xzt*r.wyt+2*r.xyzt*r.wt-2*r.wx*r.yz-2*r.w*r.xyz+2*r.wxt*r.yzt
		+2*r.z*r.wxy+2*r.xz*r.wy-2*r.wxyt*r.zt-2*r.wxz*r.y-2*r.wz*r.xy+2*r.yt*r.wxzt
		-2*r.xt*r.wyzt+2*r.q*r.wxyz-2*r.t*r.wxyzt+2*r.wzt*r.xyt+2*r.x*r.wyz ; 

	E = 2*r.wzt*r.xyz-2*r.yt*r.wx+2*r.wxzt*r.yz+2*r.t*r.wxy+2*r.wxz*r.yzt
		+2*r.xt*r.wy+2*r.wxyt*r.q-2*r.zt*r.wxyz-2*r.wxt*r.y+2*r.x*r.wyt+2*r.xyzt*r.wz
		-2*r.z*r.wxyzt-2*r.w*r.xyt-2*r.xzt*r.wyz-2*r.wt*r.xy-2*r.xz*r.wyzt ; 

	D = -2*r.w*r.xzt+2*r.wxzt*r.q-2*r.wxyt*r.yz+2*r.wxyzt*r.y-2*r.yzt*r.wxy
		+2*r.wyz*r.xyt+2*r.wyzt*r.xy-2*r.xyz*r.wyt-2*r.z*r.wxt-2*r.xz*r.wt+2*r.wxz*r.t
		+2*r.wz*r.xt-2*r.zt*r.wx+2*r.wzt*r.x-2*r.xyzt*r.wy+2*r.yt*r.wxyz ; 

	C = 2*r.t*r.wyz-2*r.wxzt*r.xy+2*r.wzt*r.y+2*r.xyzt*r.wx+2*r.wyzt*r.q
		-2*r.w*r.yzt-2*r.xt*r.wxyz+2*r.wxt*r.xyz-2*r.zt*r.wy-2*r.wxyzt*r.x+2*r.xzt*r.wxy
		-2*r.z*r.wyt+2*r.yt*r.wz-2*r.wxz*r.xyt+2*r.wxyt*r.xz-2*r.wt*r.yz ; 

	B = 2*r.yt*r.xz-2*r.wxt*r.wyz-2*r.zt*r.xy-2*r.wxyt*r.wz+2*r.t*r.xyz
		-2*r.xt*r.yz-2*r.wzt*r.wxy+2*r.wxzt*r.wy+2*r.xzt*r.y+2*r.xyzt*r.q-2*r.x*r.yzt
		-2*r.wx*r.wyzt+2*r.w*r.wxyzt+2*r.wt*r.wxyz-2*r.z*r.xyt+2*r.wxz*r.wyt ; 

	A = 2*r.xzt*r.wy+2*r.wxz*r.yt+2*r.xyzt*r.w+2*r.wt*r.xyz+2*r.xz*r.wyt
		-2*r.z*r.wxyt-2*r.yzt*r.wx-2*r.xt*r.wyz-2*r.wzt*r.xy-2*r.wxt*r.yz+2*r.t*r.wxyz
		-2*r.x*r.wyzt+2*r.wxyzt*r.q-2*r.zt*r.wxy+2*r.wxzt*r.y-2*r.wz*r.xyt ; 

	// MV2 = Reverse_Vector_Quad(r*Reverse(r));
	MV2 = Set(a, -b,-c,-d,-e,-f, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, -F,-E,-D,-C,-B, A);

	g = a*a - b*b -c*c - d*d - e*e + f*f - F*F + E*E + D*D + C*C + B*B - A*A;

	G = + 2*a*A - 2*b*B + 2*c*C - 2*d*D + 2*e*E - 2*f*F ; 

	MV3 = Set(g, 0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0, -G);

	det2 = g*g + G*G;
	
	Result = Product(Product(MV1,MV2),MV3);	// need to divide by det2. (Result is the Adjugate.)
	Result.q = Result.q/det2;

	Result.w = Result.w/det2;
	Result.x = Result.x/det2;
	Result.y = Result.y/det2;
	Result.z = Result.z/det2;
	Result.t = Result.t/det2;

	Result.wx = Result.wx/det2;
	Result.wy = Result.wy/det2;
	Result.wz = Result.wz/det2;
	Result.wt = Result.wt/det2;
	Result.xy = Result.xy/det2;
	Result.xz = Result.xz/det2;
	Result.yz = Result.yz/det2;
	Result.xt = Result.xt/det2;
	Result.yt = Result.yt/det2;
	Result.zt = Result.zt/det2;

	Result.wxy = Result.wxy/det2;
	Result.wxz = Result.wxz/det2;
	Result.wxt = Result.wxt/det2;
	Result.wyz = Result.wyz/det2;
	Result.wyt = Result.wyt/det2;
	Result.wzt = Result.wzt/det2;
	Result.xyz = Result.xyz/det2;
	Result.xyt = Result.xyt/det2;
	Result.xzt = Result.xzt/det2;
	Result.yzt = Result.yzt/det2;

	Result.wxyz = Result.wxyz/det2;
	Result.wxyt = Result.wxyt/det2;
	Result.wxzt = Result.wxzt/det2;
	Result.wyzt = Result.wyzt/det2;
	Result.xyzt = Result.xyzt/det2;

	Result.wxyzt = Result.wxyzt/det2;

	return Result;
}

//	int	main(void) { return 0; }

