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

/////////////////////////////////// Test OverBar


int main(void)
{
	
// 	ex q,  x,y,z,w, xy,xz,yz,xw,yw,zw,  xyz,xyw,xzw,yzw,  xyzw;

	GA4E q,  x,y,z,w, xy,xz,yz,xw,yw,zw,  xyz,xyw,xzw,yzw,  xyzw;

	// initialize the convenience variables above.
	
	q.q = 1;	x.x = 1;	y.y = 1;	z.z = 1;	w.w = 1;	
	xy.xy = 1;	xz.xz = 1;	yz.yz = 1;	xw.xw = 1;	yw.yw = 1;	zw.zw = 1;
	xyz.xyz = 1;	xyw.xyw = 1;	xzw.xzw = 1;	yzw.yzw = 1;	xyzw.xyzw = 1;


	GA4E Blade[16];

	GA4E MV1, MV2, MV3;

	int i, count;
	i = 0;
	count = 0;
	Blade[i++] = q; 

	Blade[i++] = x;
	Blade[i++] = y;
	Blade[i++] = z;
	Blade[i++] = w; 

	Blade[i++] = xy;
	Blade[i++] = xz;
	Blade[i++] = yz;
	Blade[i++] = xw;
	Blade[i++] = yw;
	Blade[i++] = zw;

	Blade[i++] = xyz;
	Blade[i++] = xyw;
	Blade[i++] = xzw;
	Blade[i++] = yzw;

	Blade[i++] = xyzw;



// OverBar(Blade[i]) = Reciprocal(Blade[i])*xyzw ;

// UnderBar(Blade[i]) = xyzw*Reciprocal(Blade[i]) ;

	count = 0;
	for (i=0;i<16;i++) {
		MV1 = OverBar(Blade[i]) ;
		MV2 = Reciprocal(Blade[i])*xyzw ;
		if (MV1==MV2) count++; else cout << "OverBar fail at i = " << i << "\n";
	}
	cout << "OverBar count = " << count << "\n";

	count = 0;
	for (i=0;i<16;i++) {
		MV1 = UnderBar(Blade[i]) ;
		MV2 = xyzw*Reciprocal(Blade[i]) ;
		if (MV1==MV2) count++; else cout << "UnderBar fail at i = " << i << "\n";
	}
	cout << "UnderBar count = " << count << "\n";



////////////////////////////////////

	return 0;
}


/*
	GA4E q, w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt;
	GA4E wxy, wxz, wxt, wyz, wyt, wzt, xyz, xyt, xzt, yzt;
	GA4E wxyz, wxyt, wxzt, wyzt, xyzt, wxyzt;

	// initialize the convenience variables above.
	
	q.q = 1;	w.w = 1;	x.x = 1;	y.y = 1;	z.z = 1;	t.t = 1;

	wx.wx = 1;	wy.wy = 1;	wz.wz = 1;	wt.wt = 1;	xy.xy = 1;
	xz.xz = 1;	xt.xt = 1;	yz.yz = 1;	yt.yt = 1;	zt.zt = 1;

	wxy.wxy = 1;	wxz.wxz = 1;	wxt.wxt = 1;	wyz.wyz = 1;	wyt.wyt = 1;
	wzt.wzt = 1;	xyz.xyz = 1;	xyt.xyt = 1;	xzt.xzt = 1;	yzt.yzt = 1;

	wxyz.wxyz = 1;	wxyt.wxyt = 1;	wxzt.wxzt = 1;	wyzt.wyzt = 1;	xyzt.xyzt = 1;
	wxyzt.wxyzt = 1;

	Blade[i++] = q; 
	Blade[i++] = w; 
	Blade[i++] = x;
	Blade[i++] = y;
	Blade[i++] = z;
	Blade[i++] = t;
	Blade[i++] = wx;
	Blade[i++] = wy;
	Blade[i++] = wz;
	Blade[i++] = wt;
	Blade[i++] = xy;
	Blade[i++] = xz;
	Blade[i++] = xt;
	Blade[i++] = yz;
	Blade[i++] = yt;
	Blade[i++] = zt;
	Blade[i++] = wxy;
	Blade[i++] = wxz;
	Blade[i++] = wxt;
	Blade[i++] = wyz;
	Blade[i++] = wyt;
	Blade[i++] = wzt;
	Blade[i++] = xyz;
	Blade[i++] = xyt;
	Blade[i++] = xzt;
	Blade[i++] = yzt;
	Blade[i++] = wxyz;
	Blade[i++] = wxyt;
	Blade[i++] = wxzt;
	Blade[i++] = wyzt;
	Blade[i++] = xyzt;
	Blade[i++] = wxyzt;

*/


