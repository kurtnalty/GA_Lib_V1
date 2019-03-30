#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

#include "GA5_4_1_Routines.cp"

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

int main(void)
{
	

	symbol A("A");
	symbol B("B"), C("C"), D("D"), E("E");
	symbol F("F"), G("G"), H("H"), J("J"), K("K"), L("L");
	symbol M("M"), N("N"), P("P"), R("R");
	symbol S("S");


	symbol a("a");
	symbol b("b"), c("c"), d("d"), e("e");
	symbol f("f"), g("g"), h("h"), j("j"), k("k"), l("l");
	symbol m("m"), n("n"), p("p"), r("r");
	symbol s("s");

	symbol a_q("a.q");
	symbol a_w("a.w"), a_x("a.x"), a_y("a.y"), a_z("a.z"), a_t("a.t");
	symbol a_wx("a.wx"), a_wy("a.wy"), a_wz("a.wz"), a_wt("a.wt"), a_xy("a.xy");
	symbol a_xz("a.xz"), a_xt("a.xt"), a_yz("a.yz"), a_yt("a.yt"), a_zt("a.zt");
	symbol a_wxy("a.wxy"), a_wxz("a.wxz"), a_wxt("a.wxt"), a_wyz("a.wyz"), a_wyt("a.wyt");
	symbol a_wzt("a.wzt"), a_xyz("a.xyz"), a_xyt("a.xyt"), a_xzt("a.xzt"), a_yzt("a.yzt");
	symbol a_wxyz("a.wxyz"), a_wxyt("a.wxyt"), a_wxzt("a.wxzt"), a_wyzt("a.wyzt"), a_xyzt("a.xyzt");
	symbol a_wxyzt("a.wxyzt");

	symbol b_q("b.q");
	symbol b_w("b.w"), b_x("b.x"), b_y("b.y"), b_z("b.z"), b_t("b.t");
	symbol b_wx("b.wx"), b_wy("b.wy"), b_wz("b.wz"), b_wt("b.wt"), b_xy("b.xy");
	symbol b_xz("b.xz"), b_xt("b.xt"), b_yz("b.yz"), b_yt("b.yt"), b_zt("b.zt");
	symbol b_wxy("b.wxy"), b_wxz("b.wxz"), b_wxt("b.wxt"), b_wyz("b.wyz"), b_wyt("b.wyt");
	symbol b_wzt("b.wzt"), b_xyz("b.xyz"), b_xyt("b.xyt"), b_xzt("b.xzt"), b_yzt("b.yzt");
	symbol b_wxyz("b.wxyz"), b_wxyt("b.wxyt"), b_wxzt("b.wxzt"), b_wyzt("b.wyzt"), b_xyzt("b.xyzt");
	symbol b_wxyzt("b.wxyzt");

	symbol c_q("c.q");
	symbol c_w("c.w"), c_x("c.x"), c_y("c.y"), c_z("c.z"), c_t("c.t");
	symbol c_wx("c.wx"), c_wy("c.wy"), c_wz("c.wz"), c_wt("c.wt"), c_xy("c.xy");
	symbol c_xz("c.xz"), c_xt("c.xt"), c_yz("c.yz"), c_yt("c.yt"), c_zt("c.zt");
	symbol c_wxy("c.wxy"), c_wxz("c.wxz"), c_wxt("c.wxt"), c_wyz("c.wyz"), c_wyt("c.wyt");
	symbol c_wzt("c.wzt"), c_xyz("c.xyz"), c_xyt("c.xyt"), c_xzt("c.xzt"), c_yzt("c.yzt");
	symbol c_wxyz("c.wxyz"), c_wxyt("c.wxyt"), c_wxzt("c.wxzt"), c_wyzt("c.wyzt"), c_xyzt("c.xyzt");
	symbol c_wxyzt("c.wxyzt");

	GA5_4_1 MV1,MV2,MV3,MV4,MV5;

	int i;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

//	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
//		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
//		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
//		107,109,113,127,131,   137) ; 

/*
	CA = BC   => B = C A C^{-1}    , det(B) = det(A) for these transformations

	A = MV1     B = Reverse(MV1)    C = unknown to be found

	First step, component equations for CA - BC =  0
*/
///////////////////////////////////////

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	MV2 = Reverse(MV1);
	MV3 = GA5_4_1(c_q, c_w,c_x,c_y,c_z,c_t, 
		c_wx, c_wy, c_wz, c_wt, c_xy, c_xz, c_xt, c_yz, c_yt, c_zt, 
		c_wxy, c_wxz, c_wxt, c_wyz, c_wyt, c_wzt, c_xyz, c_xyt, c_xzt, c_yzt,
		c_wxyz, c_wxyt, c_wxzt, c_wyzt, c_xyzt, c_wxyzt);

	MV4 = MV3*MV1 - MV2*MV3;

	cout << MV4.q     << " = 0 \n\n";

	cout << MV4.w     << " = 0 \n\n";
	cout << MV4.x     << " = 0 \n\n";
	cout << MV4.y     << " = 0 \n\n";
	cout << MV4.z     << " = 0 \n\n";
	cout << MV4.t     << " = 0 \n\n";


	cout << MV4.wx    << " = 0 \n\n";
	cout << MV4.wy    << " = 0 \n\n";
	cout << MV4.wz    << " = 0 \n\n";
	cout << MV4.wt    << " = 0 \n\n";
	cout << MV4.xy    << " = 0 \n\n";
	cout << MV4.xz    << " = 0 \n\n";
	cout << MV4.xt    << " = 0 \n\n";
	cout << MV4.yz    << " = 0 \n\n";
	cout << MV4.yt    << " = 0 \n\n";
	cout << MV4.zt    << " = 0 \n\n";

	cout << MV4.wxy   << " = 0 \n\n";
	cout << MV4.wxz   << " = 0 \n\n";
	cout << MV4.wxt   << " = 0 \n\n";
	cout << MV4.wyz   << " = 0 \n\n";
	cout << MV4.wyt   << " = 0 \n\n";
	cout << MV4.wzt   << " = 0 \n\n";
	cout << MV4.xyz   << " = 0 \n\n";
	cout << MV4.xyt   << " = 0 \n\n";
	cout << MV4.xzt   << " = 0 \n\n";
	cout << MV4.yzt   << " = 0 \n\n";

	cout << MV4.wxyz  << " = 0 \n\n";
	cout << MV4.wxyt  << " = 0 \n\n";
	cout << MV4.wxzt  << " = 0 \n\n";
	cout << MV4.wyzt  << " = 0 \n\n";
	cout << MV4.xyzt  << " = 0 \n\n";

	cout << MV4.wxyzt << " = 0 \n\n";

///////////////////////////////////////


	GA5_4_1 q, w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt;
	GA5_4_1 wxy, wxz, wxt, wyz, wyt, wzt, xyz, xyt, xzt, yzt;
	GA5_4_1 wxyz, wxyt, wxzt, wyzt, xyzt, wxyzt;

	// initialize the convenience variables above.
	
	q.q = 1; 
	w.w = 1; 
	x.x = 1;
	y.y = 1;
	z.z = 1;
	t.t = 1;
	wx.wx = 1;
	wy.wy = 1;
	wz.wz = 1;
	wt.wt = 1;
	xy.xy = 1;
	xz.xz = 1;
	xt.xt = 1;
	yz.yz = 1;
	yt.yt = 1;
	zt.zt = 1;
	wxy.wxy = 1;
	wxz.wxz = 1;
	wxt.wxt = 1;
	wyz.wyz = 1;
	wyt.wyt = 1;
	wzt.wzt = 1;
	xyz.xyz = 1;
	xyt.xyt = 1;
	xzt.xzt = 1;
	yzt.yzt = 1;
	wxyz.wxyz = 1;
	wxyt.wxyt = 1;
	wxzt.wxzt = 1;
	wyzt.wyzt = 1;
	xyzt.xyzt = 1;
	wxyzt.wxyzt = 1;

	GA5_4_1 Blade[32];

	i = 0;
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

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	cout << "\nTest Reciprocal of time basis for negative sign (metric) \n";
	MV2 = t;
	MV3 = Reciprocal(MV2);
	cout << MV2 << "\n";
	cout << MV3 << "\n";

	cout << "Now do Sandwich product with each blade element\n";
	cout << "Desired target is the reverse \n";
	MV2 = Reverse(MV1);
	cout << MV2 << " <- desired target \n\n";
	for(i=0;i<32;i++) {
		MV2 = Blade[i]*MV1*Reciprocal(Blade[i]);
		cout << MV2 << "\n";

		MV3 = MV1*MV2;
		if(ZCount(MV3) != 0) cout << "ZCount(MV1*MV2) = " << ZCount(MV3) << "\n";
	}

	cout << "\nTest Reciprocal of blade plus it's dual \n";
	MV2 = x + wyzt;
	MV3 = Reciprocal(MV2);
	cout << MV2 << "\n";
	cout << MV3 << "\n";
	cout << "Now check blades plus their dual \n";
	cout << "Desired target is the reverse \n";
	MV2 = Reverse(MV1);
	cout << MV2 << " <- desired target \n\n";
	for(i=0;i<16;i++) {
		MV2 = Blade[i] + Blade[31 - i];
		MV3 = Reciprocal(MV2);
		MV4 = MV2*MV1*MV3;
		cout << MV4 << "\n";

		MV5 = MV1*MV4;
		if(ZCount(MV5) != 0) cout << "ZCount(MV1*(MV2*MV1*MV3)) = " << ZCount(MV5) << "\n";
	}

	cout << "What is the determinant of idempotents and null factors? \n";
	MV2 = q/2 + x/2;
	MV3 = MV2*MV2;
	cout << "    MV2 = " << MV2 << "\n";
	cout << "MV2*MV2 = " << MV3 << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	MV4 = q/2 - x/2;
	MV5 = MV4*MV4;	// expect MV4
	cout << "    MV4 = " << MV4 << "\n";
	cout << "MV4*MV4 = " << MV5 << "\n";
	MV5 = MV2*MV4;	// expect 0
	cout << "Determinant(MV4) = " << Determinant(MV4) << "\n";
	cout << "MV2*MV4 = " << MV5 << " expect 0\n";



/*	cout << q << "\n";
	cout << w << "\n";
	cout << x << "\n";
	cout << y << "\n";
	cout << z << "\n";
	cout << t << "\n";
	cout << wx << "\n";
	cout << wy << "\n";
	cout << wz << "\n";
	cout << wt << "\n";
	cout << xy << "\n";
	cout << xz << "\n";
	cout << xt << "\n";
	cout << yz << "\n";
	cout << yt << "\n";
	cout << zt << "\n";
	cout << wxy << "\n";
	cout << wxz << "\n";
	cout << wxt << "\n";
	cout << wyz << "\n";
	cout << wyt << "\n";
	cout << wzt << "\n";
	cout << xyz << "\n";
	cout << xyt << "\n";
	cout << xzt << "\n";
	cout << yzt << "\n";
	cout << wxyz << "\n";
	cout << wxyt << "\n";
	cout << wxzt << "\n";
	cout << wyzt << "\n";
	cout << xyzt << "\n";
	cout << wxyzt << "\n";
*/

/////////////////////////////////

	return 0;
}


