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

/////////////////////////////////// Demonstrate GA5_4_1 Routines ////////////////////////////////


int main(void)
{
	
	GA5_4_1 Zero;   Zero = GA5_4_1();


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



	GA5_4_1 MV1,MV2,MV3,MV4;

	GA5_4_1 X;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);



	X = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	ex det_ref, det;
	det_ref = Determinant(MV3);

//////////////////////////////////////////////////////

// ostream &operator<<(ostream &ff, GA5_4_1 &v) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	cout << "\nMV1 = " << MV1 << "\n";

//////////////////////////////////////////////////////

// void PrintMV(GA5_4_1 &v) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	cout << "MV1 = "; PrintMV(MV1); cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 OverBar(GA5_4_1 a) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = OverBar(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "OverBar(MV1) = "; PrintMV(MV3); cout << "\n";
	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*OverBar(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 UnderBar(GA5_4_1 a) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = UnderBar(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "UnderBar(MV1) = "; PrintMV(MV3); cout << "\n";
	MV2 = OverBar(MV3);
	cout << "OverBar(UnderBar(MV1)) = "; PrintMV(MV2); cout << "\n";
	MV3 = OverBar(MV1);	MV2 = UnderBar(MV3);
	cout << "UnderBar(OverBar(MV1)) = "; PrintMV(MV2); cout << "\n";

	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*UnderBar(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 Reverse(GA5_4_1 w) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = Reverse(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "Reverse(MV1) = "; PrintMV(MV3); cout << "\n";
	MV2 = Reverse(MV3);
	cout << "Reverse(Reverse(MV1)) = "; PrintMV(MV2); cout << "\n";

	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*Reverse(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 Transpose(GA5_4_1 w) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = Transpose(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "Transpose(MV1) = "; PrintMV(MV3); cout << "\n";
	MV2 = Transpose(MV3);
	cout << "Transpose(Transpose(MV1)) = "; PrintMV(MV2); cout << "\n";
	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*Transpose(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 Conjugation(GA5_4_1 w) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = Conjugation(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "Conjugation(MV1) = "; PrintMV(MV3); cout << "\n";
	MV2 = Conjugation(MV3);
	cout << "Conjugation(Conjugation(MV1)) = "; PrintMV(MV2); cout << "\n";
	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*Conjugation(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 CliffordConjugation(GA5_4_1 w) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = CliffordConjugation(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "CliffordConjugation(MV1) = "; PrintMV(MV3); cout << "\n";
	MV2 = CliffordConjugation(MV3);
	cout << "CliffordConjugation(CliffordConjugation(MV1)) = "; PrintMV(MV2); cout << "\n";
	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*CliffordConjugation(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 Dual(GA5_4_1 w) ; 

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = Dual(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "Dual(MV1) = "; PrintMV(MV3); cout << "\n";
	MV2 = Dual(MV3);
	cout << "Dual(Dual(MV1)) = "; PrintMV(MV2); cout << "\n";
	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*Dual(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 DorstDual(GA5_4_1 a) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = DorstDual(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "DorstDual(MV1) = "; PrintMV(MV3); cout << "\n";
	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*DorstDual(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 DorstUnDual(GA5_4_1 a) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	MV3 = DorstUnDual(MV1);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "DorstUnDual(MV1) = "; PrintMV(MV3); cout << "\n";
	cout << "\n";

	MV2 = DorstDual(MV1);
	MV3 = DorstUnDual(MV2);

	cout << "DorstUnDual(DorstDual(MV1)) = "; PrintMV(MV3); cout << "\n";

	MV2 = DorstUnDual(MV1);
	MV3 = DorstDual(MV2);

	cout << "DorstDual(DorstUnDual(MV1)) = "; PrintMV(MV3); cout << "\n";

	MV3 = DorstUnDual(MV1);
	det_ref = Determinant(MV1);
	det = Determinant(MV3);
	if(det == det_ref) cout << "Determinant is conserved\n"; else cout << "Determinant is not conserved\n";
	MV2 = MV1*MV3;
	cout << "MV1*DorstUnDual(MV1) = " << MV2  << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 operator+(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = MV1 + MV2;

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV1 + MV2 = "; PrintMV(MV3); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 operator-(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = MV1 - MV2;

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV1 - MV2 = "; PrintMV(MV3); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator==(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = MV2; 

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV3 = "; PrintMV(MV3); cout << "\n";
	if (MV1 == MV2) cout << "MV1 == MV2\n";   else cout << "MV1 != MV2 (as expected)\n";
	if (MV2 == MV3) cout << "MV2 == MV3 (as expected)\n";   else cout << "MV2 != MV3 \n";
	cout << "\n";

//////////////////////////////////////////////////////

// int operator!=(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = MV2; 

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV3 = "; PrintMV(MV3); cout << "\n";
	if (MV1 != MV2) cout << "MV1 != MV2 (as expected)\n";   else cout << "MV1 == MV2 \n";
	if (MV2 != MV3) cout << "MV2 != MV3 \n";   else cout << "MV2 == MV3 (as expected)\n";
	cout << "\n";

//////////////////////////////////////////////////////

GA5_4_1 operator*(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = MV1*MV2;
	cout << "MV3 = MV1*MV2 = " << MV3 << "\n";

	MV3 = MV1*MV2 - MV2*MV1;
	if (MV3 == Zero) cout << "Product is commutative\n"; else cout << "Product is non-commutative (as expected)\n";
	cout << "\n";

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
 
	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV3 = "; PrintMV(MV3); cout << "\n";
	MV4 = (MV1*MV2)*MV3 - MV1*(MV2*MV3);
	if (MV4 == Zero) cout << "Product is associative (as expected)\n"; else cout << "Product is non-associative\n";
	cout << "\n";

	MV4 = symbol("Z")*MV3;
	cout << "MV4 = Z*MV3 = "; PrintMV(MV4); cout << "\n";
	MV4 = MV3*symbol("Z");
	cout << "MV4 = MV3*Z = "; PrintMV(MV4); cout << "\n";
	cout << "\n";


//////////////////////////////////////////////////////

// GA5_4_1 Product(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Product(MV1,MV2);
	cout << "MV3 = Product(MV1,MV2) = " << MV3 << "\n";

	MV3 = Product(MV1,MV2) - Product(MV2,MV1);
	if (MV3 == Zero) cout << "Product(MV1,MV2) is commutative\n"; else cout << "Product(MV1,MV2) is non-commutative (as expected)\n";
	cout << "\n";

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
 
	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV3 = "; PrintMV(MV3); cout << "\n";
	MV4 = Product(Product(MV1,MV2), MV3) - Product(MV1,Product(MV2,MV3));
	if (MV4 == Zero) cout << "Product is associative (as expected)\n"; else cout << "Product is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 operator/(const GA5_4_1 &u, const int i) ;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);  
	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	MV2 = MV1/2;
	cout << "MV2 = MV1/2 = "; PrintMV(MV2); cout << "\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 operator^(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = MV1^MV2;
	cout << "MV3 = MV1^MV2 = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 Wedge(const GA5_4_1 &u, const GA5_4_1 &v) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Wedge(MV1,MV2);
	cout << "MV3 = Wedge(MV1,MV2) = " << MV3 << "\n";

	MV3 = Wedge(MV1,MV2) - Wedge(MV2,MV1);
	if (MV3 == Zero) cout << "Wedge(MV1,MV2) is commutative\n"; else cout << "Wedge(MV1,MV2) is non-commutative (as expected)\n";
	cout << "\n";

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A); 
	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV3 = "; PrintMV(MV3); cout << "\n";
	MV4 = Wedge(Wedge(MV1,MV2), MV3) - Wedge(MV1,Wedge(MV2,MV3));
	if (MV4 == Zero) cout << "Wedge is associative (as expected)\n"; else cout << "Wedge is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 AntiWedge(const GA5_4_1 a, const GA5_4_1 b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = AntiWedge(MV1,MV2);
	cout << "MV3 = AntiWedge(MV1,MV2) = " << MV3 << "\n";

	MV3 = AntiWedge(MV1,MV2) - AntiWedge(MV2,MV1);
	if (MV3 == Zero) cout << "AntiWedge(MV1,MV2) is commutative\n"; else cout << "AntiWedge(MV1,MV2) is non-commutative (as expected)\n";
	cout << "\n";

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A); 
	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV3 = "; PrintMV(MV3); cout << "\n";
	MV4 = AntiWedge(AntiWedge(MV1,MV2), MV3) - AntiWedge(MV1,AntiWedge(MV2,MV3));
	if (MV4 == Zero) cout << "AntiWedge is associative (as expected)\n"; else cout << "AntiWedge is non-associative\n";
	cout << "\n";


//////////////////////////////////////////////////////

// GA5_4_1 Regressive(GA5_4_1 a, GA5_4_1 b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Regressive(MV1,MV2);
	cout << "MV3 = Regressive(MV1,MV2) = " << MV3 << "\n";

	MV3 = Regressive(MV1,MV2) - Regressive(MV2,MV1);
	if (MV3 == Zero) cout << "Regressive(MV1,MV2) is commutative\n"; else cout << "Regressive(MV1,MV2) is non-commutative (as expected)\n";
	cout << "\n";

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV3 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A); 
	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	cout << "MV3 = "; PrintMV(MV3); cout << "\n";
	MV4 = Regressive(Regressive(MV1,MV2), MV3) - Regressive(MV1,Regressive(MV2,MV3));
	if (MV4 == Zero) cout << "Regressive is associative (as expected)\n"; else cout << "Regressive is non-associative\n";
	cout << "\n";

//////////////////////////////////////////////////////

// GA5_4_1 RegressiveViaFormula(GA5_4_1 a, GA5_4_1 b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = RegressiveViaFormula(MV1,MV2);
	cout << "MV3 = RegressiveViaFormula(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 LowerRightViaFormula(GA5_4_1 a, GA5_4_1 b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = LowerRightViaFormula(MV1,MV2);
	cout << "MV3 = LowerRightViaFormula(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 Expander(const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Expander(MV1,MV2);
	cout << "MV3 = Expander(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 Conserver(const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Conserver(MV1,MV2);
	cout << "MV3 = Conserver(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 Shrinker(const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Shrinker(MV1,MV2);
	cout << "MV3 = Shrinker(MV1,MV2) = " << MV3 << "\n\n";

///////////////////////////////////////////////////////

// GA5_4_1 Symmetric(const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Symmetric(MV1,MV2);
	cout << "MV3 = Symmetric(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 AntiSymmetric(const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = AntiSymmetric(MV1,MV2);
	cout << "MV3 = AntiSymmetric(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 Inner(const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = Inner(MV1,MV2);
	cout << "MV3 = Inner(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 LeftContraction (const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = LeftContraction(MV1,MV2);
	cout << "MV3 = LeftContraction(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// GA5_4_1 RightContraction (const GA5_4_1 &a, const GA5_4_1 &b) ;

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	cout << "MV1 = "; PrintMV(MV1); cout << "\n";
	cout << "MV2 = "; PrintMV(MV2); cout << "\n";
	MV3 = RightContraction(MV1,MV2);
	cout << "MV3 = RightContraction(MV1,MV2) = " << MV3 << "\n\n";

//////////////////////////////////////////////////////

// ex Determinant(GA5_4_1 A) ;

	X = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	cout << "X = "; PrintMV(X); cout << "\n";
	cout << "Determinant(X) = " << Determinant(X) << "\n\n";


//////////////////////////////////////////////////////

// GA5_4_1 Adjugate(GA5_4_1 V) ;

	X = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	cout << "X = "; PrintMV(X); cout << "\n";
	MV3 = Adjugate(X);
	cout << "MV3 = Adjugate(X) = " << MV3 << "\n\n";
	MV4 = MV3*X;
	cout << "Adjugate(X)*X = " << MV4 << "\n\n";
	MV4 = X*MV3;
	cout << "X*Adjugate(X) = " << MV4 << "\n\n";


//////////////////////////////////////////////////////

// GA5_4_1 Reciprocal(GA5_4_1 a) ;

	X = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	cout << "X = "; PrintMV(X); cout << "\n";
	MV3 = Reciprocal(X);
	cout << "MV3 = Reciprocal(X) = " << MV3 << "\n";
	MV2 = X*MV3;
	cout << "X*(1/X) = " << MV2 << "\n\n";

/////////////////////////////////  

// now I want to look at the complex determinant for a prime number test vector for the mismatched unary operations

	cout << "\nUnary Operator Printout\n\n";

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	cout << "Reference input MV for each function below \n" << MV1 << "\n\n";

	MV2 = OverBar(MV1) ;  
	cout << "OverBar \n" << MV2 << "\n\n";

	MV2 = UnderBar(MV1) ;  
	cout << "UnderBar \n" << MV2 << "\n\n";

	MV2 = Reverse(MV1) ;  
	cout << "Reverse \n" << MV2 << "\n\n";

//	MV2 = Reverse_Vector_Quad(MV1) ;  
//	cout << "Reverse_Vector_Quad \n" << MV2 << "\n\n";

	MV2 = Transpose(MV1) ;  
	cout << "Transpose \n" << MV2 << "\n\n";

//	MV2 = Conjugation(MV1) ;  
//	cout << "Conjugation \n" << MV2 << "\n\n";

	MV2 = CliffordConjugation(MV1) ;  
	cout << "CliffordConjugation \n" << MV2 << "\n\n";

	MV2 = Dual(MV1) ;  
	cout << "Dual \n" << MV2 << "\n\n";

	MV2 = DorstDual(MV1) ;  
	cout << "DorstDual \n" << MV2 << "\n\n";

	MV2 = DorstUnDual(MV1) ;  
	cout << "DorstUnDual \n" << MV2 << "\n\n";

	MV2 = Parity(MV1) ;  
	cout << "Parity \n" << MV2 << "\n\n";

	MV2 = ComplexConjugate(MV1) ;  
	cout << "ComplexConjugate \n" << MV2 << "\n\n";

	MV2 = Hermitian(MV1) ;  
	cout << "Hermitian \n" << MV2 << "\n\n";

//	MV2 = Conjugate(MV1) ;  
//	cout << "Conjugate \n" << MV2 << "\n\n";
	 

	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	cout << "\nNotice the complex conjugates below. For pure magnitudes, these determinants would agree.\n" ;
	cout << "I am keeping determinant complex to allow easy comparison to Gamma matrix results in QM.\n\n" ;

	cout << "Reference \n" << Determinant(MV1) << "\n\n";

	MV2 = Reverse(MV1) ;  
	cout << "Reverse \n" << Determinant(MV2) << "\n\n";

	MV2 = Transpose(MV1) ;  
	cout << "Transpose \n" << Determinant(MV2) << "\n\n";

	MV2 = Dual(MV1) ;  
	cout << "Dual \n" << Determinant(MV2) << "\n\n";

	MV2 = DorstDual(MV1) ;  
	cout << "DorstDual \n" << Determinant(MV2) << "\n\n";

	MV2 = DorstUnDual(MV1) ;  
	cout << "DorstUnDual \n" << Determinant(MV2) << "\n\n";
///////////////////////
	MV2 = OverBar(MV1) ;  
	cout << "OverBar \n" << Determinant(MV2) << "\n\n";

	MV2 = UnderBar(MV1) ;  
	cout << "UnderBar \n" << Determinant(MV2) << "\n\n";

//	MV2 = Reverse_Vector_Quad(MV1) ;  
//	cout << "Reverse_Vector_Quad \n" << Determinant(MV2) << "\n\n";

//	MV2 = Conjugation(MV1) ;  
//	cout << "Conjugation \n" << Determinant(MV2) << "\n\n";

	MV2 = CliffordConjugation(MV1) ;  
	cout << "CliffordConjugation \n" << Determinant(MV2) << "\n\n";

	MV2 = Parity(MV1) ;  
	cout << "Parity \n" << Determinant(MV2) << "\n\n";

	MV2 = ComplexConjugate(MV1) ;  
	cout << "ComplexConjugate \n" << Determinant(MV2) << "\n\n";

	MV2 = Hermitian(MV1) ;  
	cout << "Hermitian \n" << Determinant(MV2) << "\n\n";

//	MV2 = Conjugate(MV1) ;  
//	cout << "Conjugate \n" << Determinant(MV2) << "\n\n";
	 
	cout << "\nNow we show various products of a MV and the named unary operator.\n" ;
	cout << "The reverse product is the clear winner in terms of zeroed product terms, and faithfully reserves the determinant.\n\n" ;

	MV3 = MV1*OverBar(MV1) ;
	cout << "MV1*OverBar(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*UnderBar(MV1) ;
	cout << "MV1*UnderBar(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*CliffordConjugation(MV1) ;
	cout << "MV1*CliffordConjugation(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*Parity(MV1) ;
	cout << "MV1*Parity(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*ComplexConjugate(MV1) ;
	cout << "MV1*ComplexConjugate(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*Hermitian(MV1) ;
	cout << "MV1*Hermitian(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*Reverse(MV1) ;
	cout << "MV1*Reverse(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*Transpose(MV1) ;
	cout << "MV1*Transpose(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*Dual(MV1) ;
	cout << "MV1*Dual(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*DorstUnDual(MV1) ;
	cout << "MV1*DorstUnDual(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = Reverse(MV1)*Hermitian(MV1) ;
	cout << "Reverse(MV1)*Hermitian(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

//////////////////////////////////////////

	cout << "\n\nCheck commutativity of the pseudo-scalar wxyzt \n\n";



	GA5_4_1 q, w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt;
	GA5_4_1 wxy, wxz, wxt, wyz, wyt, wzt, xyz, xyt, xzt, yzt;
	GA5_4_1 wxyz, wxyt, wxzt, wyzt, xyzt, wxyzt;

	// initialize the convenience variables above.
	
	q.q = 1;	w.w = 1;	x.x = 1;	y.y = 1;	z.z = 1;	t.t = 1;

	wx.wx = 1;	wy.wy = 1;	wz.wz = 1;	wt.wt = 1;	xy.xy = 1;
	xz.xz = 1;	xt.xt = 1;	yz.yz = 1;	yt.yt = 1;	zt.zt = 1;

	wxy.wxy = 1;	wxz.wxz = 1;	wxt.wxt = 1;	wyz.wyz = 1;	wyt.wyt = 1;
	wzt.wzt = 1;	xyz.xyz = 1;	xyt.xyt = 1;	xzt.xzt = 1;	yzt.yzt = 1;

	wxyz.wxyz = 1;	wxyt.wxyt = 1;	wxzt.wxzt = 1;	wyzt.wyzt = 1;	xyzt.xyzt = 1;
	wxyzt.wxyzt = 1;
	symbol aa("aa"), AA("AA");	// used in complex number format for determinant, etc

	GA5_4_1 Blade[32];

	GA5_4_1 MV5;

	int i, count;
	i = 0;
	count = 0;
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

	count = 0;
	for (i=0;i<32;i++) {
		MV1 = Blade[i]*wxyzt - wxyzt*Blade[i];
		if (MV1==Zero) count++; else cout << "Commute fail at i = " << i << "\n";
	}
	cout << "Commutator count = " << count << "\n";

	count = 0;
	for (i=0;i<32;i++) {
		if (Determinant(Blade[i])==1) count++; else cout << "Non-unit blade determinant at i = " << i << "\n";
	}
	cout << "Unit Determinant Blade count = " << count << "\n";

////////////////////////////////////

// check commutativity of MV and complex number

	MV1 = GA5_4_1(a,  b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A );
	MV2 = GA5_4_1(aa, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, AA);  // complex number

	MV3 = MV1*MV2 - MV2*MV1;

	cout << "MV1 = " << MV1 << "\n";
	cout << "MV2 = " << MV2 << "\n";
	cout << "MV3 = MV1*MV2 - MV2*MV1 = " << MV3 << "\n";

////////////////////////////////////

	MV1 = GA5_4_1(a,  b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A );
	MV2 = Reverse(MV1);
	MV3 = Reverse_Vector_Quad(MV1*Reverse(MV1));
//	Adjugate = MV2*MV3;
	MV4 = MV2*MV3;
	MV5 = MV4 - Adjugate(MV1);
	cout << "Adjugate formula check - delta = " << MV5 << " expect 0 \n";


////////////////////////////////////

//	verify the 64 determinant preserving Comp transforms

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	det_ref = Determinant(MV1);

	printf("\nVerifying 64 Comp transforms\n");
	for (i=0;i<64;i++) 
{
		MV2 = Magic(MV1,i);
		det = Determinant(MV2);
		if(det != det_ref) printf("Failure at index %d \n", i);
		else printf("+ ");
}
	printf("\n");
	


////////////////////////////////////

//	verify the 360 determinant preserving Magic transforms

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	det_ref = Determinant(MV1);

	printf("\nVerifying 360 Magic transforms\n");
	for (i=0;i<360;i++) 
{
		MV2 = Magic(MV1,i);
		det = Determinant(MV2);
		if(det != det_ref) printf("Failure at index %d \n", i);
		else printf("+ ");
}
	printf("\n");
	


////////////////////////////////////

	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	matrix NM;
	NM = GA5_4_1_To_Matrix(MV1);
	MV2 = Matrix_To_GA5_4_1(NM);

	printf("\nMV1 = ");	PrintMV(MV1); 	printf("\n");
	cout << "NM = " << NM << "\n";
	printf("\nMV2 = ");	PrintMV(MV2); 	printf("\n");

////////////////////////////////////

	return 0;
}

/* copy and paste bin

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	X = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 



*/





