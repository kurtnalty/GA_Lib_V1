// Speed test search for Geometric Algebra in Five Dimensional metric(4,1) space
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
#include <time.h>

//////////////////////////////////////////////////////

typedef struct
{
	int q ;
	int w,x,y,z,t ;
	int wx,wy,wz,wt,xy,xz,xt,yz,yt,zt ;
	int wxy,wxz,wxt,wyz,wyt,wzt,xyz,xyt,xzt,yzt ;
	int wxyz,wxyt,wxzt,wyzt,xyzt; 
	int wxyzt;
} GA5_4_1;


//////////////////////////////////////////////////////


GA5_4_1 Set(int q, 
		int w, int x, int y, int z, int t, 
		int wx, int wy, int wz, int wt, int xy, 
		int xz, int xt, int yz, int yt, int zt, 
		int wxy, int wxz, int wxt, int wyz, int wyt, 
		int wzt, int xyz, int xyt, int xzt, int yzt, 
		int wxyz, int wxyt, int wxzt, int wyzt, int xyzt, 
		int wxyzt) 
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

GA5_4_1 Comp(GA5_4_1 r, long index)
{
	GA5_4_1 s;
	long Mask = 1;

	if((index & Mask) == 0) s.q     = r.q    ; else s.q     = -r.q    ;	Mask <<= 1;

	if((index & Mask) == 0) s.w     = r.w    ; else s.w     = -r.w    ;	Mask <<= 1;
	if((index & Mask) == 0) s.x     = r.x    ; else s.x     = -r.x    ;	Mask <<= 1;
	if((index & Mask) == 0) s.y     = r.y    ; else s.y     = -r.y    ;	Mask <<= 1;
	if((index & Mask) == 0) s.z     = r.z    ; else s.z     = -r.z    ;	Mask <<= 1;
	if((index & Mask) == 0) s.t     = r.t    ; else s.t     = -r.t    ;	Mask <<= 1;

	if((index & Mask) == 0) s.wx    = r.wx   ; else s.wx    = -r.wx   ;	Mask <<= 1;
	if((index & Mask) == 0) s.wy    = r.wy   ; else s.wy    = -r.wy   ;	Mask <<= 1;
	if((index & Mask) == 0) s.wz    = r.wz   ; else s.wz    = -r.wz   ;	Mask <<= 1;
	if((index & Mask) == 0) s.wt    = r.wt   ; else s.wt    = -r.wt   ;	Mask <<= 1;
	if((index & Mask) == 0) s.xy    = r.xy   ; else s.xy    = -r.xy   ;	Mask <<= 1;
	if((index & Mask) == 0) s.xz    = r.xz   ; else s.xz    = -r.xz   ;	Mask <<= 1;
	if((index & Mask) == 0) s.xt    = r.xt   ; else s.xt    = -r.xt   ;	Mask <<= 1;
	if((index & Mask) == 0) s.yz    = r.yz   ; else s.yz    = -r.yz   ;	Mask <<= 1;
	if((index & Mask) == 0) s.yt    = r.yt   ; else s.yt    = -r.yt   ;	Mask <<= 1;
	if((index & Mask) == 0) s.zt    = r.zt   ; else s.zt    = -r.zt   ;	Mask <<= 1;

	if((index & Mask) == 0) s.wxy   = r.wxy  ; else s.wxy   = -r.wxy  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxz   = r.wxz  ; else s.wxz   = -r.wxz  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxt   = r.wxt  ; else s.wxt   = -r.wxt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wyz   = r.wyz  ; else s.wyz   = -r.wyz  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wyt   = r.wyt  ; else s.wyt   = -r.wyt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.wzt   = r.wzt  ; else s.wzt   = -r.wzt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.xyz   = r.xyz  ; else s.xyz   = -r.xyz  ;	Mask <<= 1;
	if((index & Mask) == 0) s.xyt   = r.xyt  ; else s.xyt   = -r.xyt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.xzt   = r.xzt  ; else s.xzt   = -r.xzt  ;	Mask <<= 1;
	if((index & Mask) == 0) s.yzt   = r.yzt  ; else s.yzt   = -r.yzt  ;	Mask <<= 1;

	if((index & Mask) == 0) s.wxyz  = r.wxyz ; else s.wxyz  = -r.wxyz ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxyt  = r.wxyt ; else s.wxyt  = -r.wxyt ;	Mask <<= 1;
	if((index & Mask) == 0) s.wxzt  = r.wxzt ; else s.wxzt  = -r.wxzt ;	Mask <<= 1;
	if((index & Mask) == 0) s.wyzt  = r.wyzt ; else s.wyzt  = -r.wyzt ;	Mask <<= 1;
	if((index & Mask) == 0) s.xyzt  = r.xyzt ; else s.xyzt  = -r.xyzt ;	Mask <<= 1;

	if((index & Mask) == 0) s.wxyzt = r.wxyzt; else s.wxyzt = -r.wxyzt;	Mask <<= 1;

	return s;
}

//////////////////////////////////////////////////////

int	ZCount(GA5_4_1 a)
{
	int	count = 0;

	if(a.q     == 0) count++;

	if(a.w     == 0) count++;
	if(a.x     == 0) count++;
	if(a.y     == 0) count++;
	if(a.z     == 0) count++;
	if(a.t     == 0) count++;

	if(a.wx    == 0) count++;
	if(a.wy    == 0) count++;
	if(a.wz    == 0) count++;
	if(a.wt    == 0) count++;
	if(a.xy    == 0) count++;
	if(a.xz    == 0) count++;
	if(a.xt    == 0) count++;
	if(a.yz    == 0) count++;
	if(a.yt    == 0) count++;
	if(a.zt    == 0) count++;

	if(a.wxy   == 0) count++;
	if(a.wxz   == 0) count++;
	if(a.wxt   == 0) count++;
	if(a.wyz   == 0) count++;
	if(a.wyt   == 0) count++;
	if(a.wzt   == 0) count++;
	if(a.xyz   == 0) count++;
	if(a.xyt   == 0) count++;
	if(a.xzt   == 0) count++;
	if(a.yzt   == 0) count++;

	if(a.wxyz  == 0) count++;
	if(a.wxyt  == 0) count++;
	if(a.wxzt  == 0) count++;
	if(a.wyzt  == 0) count++;
	if(a.xyzt  == 0) count++;

	if(a.wxyzt == 0) count++;

	return count;
}

//////////////////////////////////////////////////////

void PrintMV(GA5_4_1 v) {

	printf("(%10d, ", v.q) ;
	printf("%10d,%10d,%10d,%10d,%10d, ", v.w, v.x, v.y, v.z, v.t) ;
	printf("%10d,%10d,%10d,%10d,%10d,%10d,%10d,%10d,%10d,%10d, ", 
		v.wx, v.wy, v.wz, v.wt, v.xy, v.xz, v.xt, v.yz, v.yt, v.zt) ;
	printf("%10d,%10d,%10d,%10d,%10d,%10d,%10d,%10d,%10d,%10d, ", 
		v.wxy, v.wxz, v.wxt, v.wyz, v.wyt, v.wzt, v.xyz, v.xyt, v.xzt, v.yzt) ;
	printf("%10d,%10d,%10d,%10d,%10d, ", v.wxyz, v.wxyt, v.wxzt, v.wyzt, v.xyzt) ;
	printf("%10d) \n", v.wxyzt);

}

//////////////////////////////////////////////////////

int main(void) 
{

	GA5_4_1 a,b,c;
	long i;

	time_t start_time, end_time;
	long deltat;

	int count;
	
	a = Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	b = Reverse(a);
	c = Product(a,b);
	count = ZCount(c);
//	PrintMV(a);
//	PrintMV(b);
//	PrintMV(c);
//	printf("count = %d \n\n",count);

//	time(&start_time);
//	for(i=0;i<1000000;i++ ) det = Determinant(r);
//	time(&end_time);
//	deltat = difftime(end_time, start_time);
//	printf("deltat for 1000000X MV determinant = %ld \n", deltat);

	time(&start_time);
	for(i=0;i < 0x80000000 ;i++ ) // Product(a,Reverse(a));
//	for(i=0;i < 0x20 ;i++ ) // Product(a,Reverse(a));
	{
		b = Comp(a,i);
		c = Product(a,b);
		count = ZCount(c);
//		PrintMV(a);
//		PrintMV(b);
//		PrintMV(c);
//		printf("count = %d \n\n",count);
		if (count > 19) printf("count = %d at index %ld \n",count, i);
	}
	time(&end_time);
	deltat = difftime(end_time, start_time);
	printf("deltat for 0x80000000 = %ld \n", deltat);

/*

deltat for 1000000X MV determinant = 1 (long)
deltat for 10000000X Product(r,Reverse(r)) = 35 


deltat for 1 000 000X MV determinant = 1 (0.67 us/determinant using greater 100X loop counter) double
deltat for 10 000 000X Product(r,Reverse(r)) = 33    3.3us/Product

deltat for 1000000X MV determinant = 0 (int)
deltat for 10000000X Product(r,Reverse(r)) = 16 converted to int from long
deltat for 10000000X Product(r,Reverse(r)) = 15 inlined Product, no subroutine calls in loop.

deltat for 10000000X Product(r,Reverse(r)) = 19 reverted to subroutines, added Comp(a,i)

deltat for 0x01000000 = 32 
*/

	return 0;
}
