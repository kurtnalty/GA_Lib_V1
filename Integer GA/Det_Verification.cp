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

	long All_Indices[128] = {
		0x00000000, 0x03FFFFC0, 0x044B2DDF, 0x07B4D21F, 0x089556EF, 0x0B6AA92F, 0x0CDE7B30, 0x0F2184F0, 
		0x11269B77, 0x12D964B7, 0x156DB6A8, 0x16924968, 0x19B3CD98, 0x1A4C3258, 0x1DF8E047, 0x1E071F87, 
		0x21C71C7B, 0x2238E3BB, 0x258C31A4, 0x2673CE64, 0x29524A94, 0x2AADB554, 0x2D19674B, 0x2EE6988B, 
		0x30E1870C, 0x331E78CC, 0x34AAAAD3, 0x37555513, 0x3874D1E3, 0x3B8B2E23, 0x3C3FFC3C, 0x3FC003FC, 
		0x403F03FD, 0x43C0FC3D, 0x44742E22, 0x478BD1E2, 0x48AA5512, 0x4B55AAD2, 0x4CE178CD, 0x4F1E870D, 
		0x5119988A, 0x52E6674A, 0x5552B555, 0x56AD4A95, 0x598CCE65, 0x5A7331A5, 0x5DC7E3BA, 0x5E381C7A, 
		0x61F81F86, 0x6207E046, 0x65B33259, 0x664CCD99, 0x696D4969, 0x6A92B6A9, 0x6D2664B6, 0x6ED99B76, 
		0x70DE84F1, 0x73217B31, 0x7495A92E, 0x776A56EE, 0x784BD21E, 0x7BB42DDE, 0x7C00FFC1, 0x7FFF0001, 
		
		0x8000FFFE, 0x83FF003E, 0x844BD221, 0x87B42DE1, 0x8895A911, 0x8B6A56D1, 0x8CDE84CE, 0x8F217B0E, 
		0x91266489, 0x92D99B49, 0x956D4956, 0x9692B696, 0x99B33266, 0x9A4CCDA6, 0x9DF81FB9, 0x9E07E079, 
		0xA1C7E385, 0xA2381C45, 0xA58CCE5A, 0xA673319A, 0xA952B56A, 0xAAAD4AAA, 0xAD1998B5, 0xAEE66775, 
		0xB0E178F2, 0xB31E8732, 0xB4AA552D, 0xB755AAED, 0xB8742E1D, 0xBB8BD1DD, 0xBC3F03C2, 0xBFC0FC02, 
		0xC03FFC03, 0xC3C003C3, 0xC474D1DC, 0xC78B2E1C, 0xC8AAAAEC, 0xCB55552C, 0xCCE18733, 0xCF1E78F3, 
		0xD1196774, 0xD2E698B4, 0xD5524AAB, 0xD6ADB56B, 0xD98C319B, 0xDA73CE5B, 0xDDC71C44, 0xDE38E384, 
		0xE1F8E078, 0xE2071FB8, 0xE5B3CDA7, 0xE64C3267, 0xE96DB697, 0xEA924957, 0xED269B48, 0xEED96488, 
		0xF0DE7B0F, 0xF32184CF, 0xF49556D0, 0xF76AA910, 0xF84B2DE0, 0xFBB4D220, 0xFC00003F, 0xFFFFFFFF
	} ;

	long Only_Matches[64] = {
		0x00000000, 0x03FFFFC0, 0x0CDE7B30, 0x0F2184F0, 0x156DB6A8, 0x16924968, 0x19B3CD98, 0x1A4C3258, 
		0x258C31A4, 0x2673CE64, 0x29524A94, 0x2AADB554, 0x30E1870C, 0x331E78CC, 0x3C3FFC3C, 0x3FC003FC, 
		0x44742E22, 0x478BD1E2, 0x48AA5512, 0x4B55AAD2, 0x5119988A, 0x52E6674A, 0x5DC7E3BA, 0x5E381C7A, 
		0x61F81F86, 0x6207E046, 0x6D2664B6, 0x6ED99B76, 0x7495A92E, 0x776A56EE, 0x784BD21E, 0x7BB42DDE, 
		0x844BD221, 0x87B42DE1, 0x8895A911, 0x8B6A56D1, 0x91266489, 0x92D99B49, 0x9DF81FB9, 0x9E07E079, 
		0xA1C7E385, 0xA2381C45, 0xAD1998B5, 0xAEE66775, 0xB4AA552D, 0xB755AAED, 0xB8742E1D, 0xBB8BD1DD, 
		0xC03FFC03, 0xC3C003C3, 0xCCE18733, 0xCF1E78F3, 0xD5524AAB, 0xD6ADB56B, 0xD98C319B, 0xDA73CE5B, 
		0xE5B3CDA7, 0xE64C3267, 0xE96DB697, 0xEA924957, 0xF0DE7B0F, 0xF32184CF, 0xFC00003F, 0xFFFFFFFF
	} ;

	long Only_Conjugates[64] = {
		0x044B2DDF, 0x07B4D21F, 0x089556EF, 0x0B6AA92F, 0x11269B77, 0x12D964B7, 0x1DF8E047, 0x1E071F87, 
		0x21C71C7B, 0x2238E3BB, 0x2D19674B, 0x2EE6988B, 0x34AAAAD3, 0x37555513, 0x3874D1E3, 0x3B8B2E23, 
		0x403F03FD, 0x43C0FC3D, 0x4CE178CD, 0x4F1E870D, 0x5552B555, 0x56AD4A95, 0x598CCE65, 0x5A7331A5, 
		0x65B33259, 0x664CCD99, 0x696D4969, 0x6A92B6A9, 0x70DE84F1, 0x73217B31, 0x7C00FFC1, 0x7FFF0001,
		0x8000FFFE, 0x83FF003E, 0x8CDE84CE, 0x8F217B0E, 0x956D4956, 0x9692B696, 0x99B33266, 0x9A4CCDA6, 
		0xA58CCE5A, 0xA673319A, 0xA952B56A, 0xAAAD4AAA, 0xB0E178F2, 0xB31E8732, 0xBC3F03C2, 0xBFC0FC02, 
		0xC474D1DC, 0xC78B2E1C, 0xC8AAAAEC, 0xCB55552C, 0xD1196774, 0xD2E698B4, 0xDDC71C44, 0xDE38E384, 
		0xE1F8E078, 0xE2071FB8, 0xED269B48, 0xEED96488, 0xF49556D0, 0xF76AA910, 0xF84B2DE0, 0xFBB4D220
	} ;

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


	GA5_4_1 MV1,MV2;
	ex det, det_ref;
	int i;

//	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 


	det_ref = Determinant(MV1);

	;
	for(i=0;i<128;i++) {
		MV2 = Comp(MV1,  All_Indices[i]);
		det = Determinant(MV2);
		if(det != det_ref) printf("0x%8lX, ", All_Indices[i]);
	}

	cout << "\n";

	return 0;
}








