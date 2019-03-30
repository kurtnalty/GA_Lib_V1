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


	GA5_4_1 MV1,MV2,MV3,MV4,MV5;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

//	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
//		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
//		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
//		107,109,113,127,131,   137) ; 

// count = 20 at index 67108800 
	MV2 = Comp(MV1,  67108800);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n\n";

// count = 20 at index 1010826300 
	MV2 = Comp(MV1,  1010826300);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n\n";

// count = 20 at index 1573381050 
	MV2 = Comp(MV1,  1573381050);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n\n";

// count = 20 at index 1859754870 
	MV2 = Comp(MV1,  1859754870);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n\n";

// count = 20 at index 2003457774 
	MV2 = Comp(MV1,  2003457774);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n\n";

// count = 20 at index 2075405790 
	MV2 = Comp(MV1,  2075405790);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n\n";

	cout << "\n\nNumerical example of the same . . \n\n";

	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

// count = 20 at index 67108800 
	MV2 = Comp(MV1,  67108800);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";

// count = 20 at index 1010826300 
	MV2 = Comp(MV1,  1010826300);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";

// count = 20 at index 1573381050 
	MV2 = Comp(MV1,  1573381050);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";

// count = 20 at index 1859754870 
	MV2 = Comp(MV1,  1859754870);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";

// count = 20 at index 2003457774 
	MV2 = Comp(MV1,  2003457774);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";

// count = 20 at index 2075405790 
	MV2 = Comp(MV1,  2075405790);
	MV3 = MV1*MV2;
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";

////////////////////////////////

	cout << "Further investigations for unit dets and ratios \n\n";

// count = 20 at index 67108800 
	MV2 = Comp(MV1,  67108800);
	MV3 = MV1*Reciprocal(MV2);
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";
	MV4 = MV2*MV3;
	MV5 = MV3*MV2;
	cout << "MV4 = " << MV4 << "\n";
	cout << "MV5 = " << MV5 << "\n";

// count = 20 at index 1010826300 
	MV2 = Comp(MV1,  1010826300);
	MV3 = MV1*Reciprocal(MV2);
	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = " << MV3 << "\n";
	cout << "Determinant(MV1) = " << Determinant(MV1) << "\n";
	cout << "Determinant(MV2) = " << Determinant(MV2) << "\n";
	cout << "Determinant(MV3) = " << Determinant(MV3) << "\n\n";


	return 0;
}


