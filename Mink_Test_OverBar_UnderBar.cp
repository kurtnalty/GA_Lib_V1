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

/////////////////////////////////// Test OverBar


int main(void)
{
	Mink q,  x,y,z,t, xy,xz,yz,xt,yt,zt,  xyz,xyt,xzt,yzt,  xyzt;

	// initialize the convenience variables above.
	
	q.q = 1;	x.x = 1;	y.y = 1;	z.z = 1;	t.t = 1;	
	xy.xy = 1;	xz.xz = 1;	yz.yz = 1;	xt.xt = 1;	yt.yt = 1;	zt.zt = 1;
	xyz.xyz = 1;	xyt.xyt = 1;	xzt.xzt = 1;	yzt.yzt = 1;	xyzt.xyzt = 1;


	Mink Blade[16];

	Mink MV1, MV2, MV3;

	int i, count;
	i = 0;
	count = 0;
	Blade[i++] = q; 

	Blade[i++] = x;
	Blade[i++] = y;
	Blade[i++] = z;
	Blade[i++] = t; 

	Blade[i++] = xy;
	Blade[i++] = xz;
	Blade[i++] = yz;
	Blade[i++] = xt;
	Blade[i++] = yt;
	Blade[i++] = zt;

	Blade[i++] = xyz;
	Blade[i++] = xyt;
	Blade[i++] = xzt;
	Blade[i++] = yzt;

	Blade[i++] = xyzt;



// OverBar(Blade[i]) = Reciprocal(Blade[i])*xyzt ;

// UnderBar(Blade[i]) = xyzt*Reciprocal(Blade[i]) ;

	count = 0;
	for (i=0;i<16;i++) {
		MV1 = OverBar(Blade[i]) ;
		MV2 = Reciprocal(Blade[i])*xyzt ;
		if (MV1==MV2) count++; else cout << "OverBar fail at i = " << i << "\n";
	}
	cout << "OverBar count = " << count << "\n";

	count = 0;
	for (i=0;i<16;i++) {
		MV1 = UnderBar(Blade[i]) ;
		MV2 = xyzt*Reciprocal(Blade[i]) ;
		if (MV1==MV2) count++; else cout << "UnderBar fail at i = " << i << "\n";
	}
	cout << "UnderBar count = " << count << "\n";



////////////////////////////////////

	return 0;
}


