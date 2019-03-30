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

///////////////////////////////////

	matrix zero(4,4);
	matrix Q(4,4);	// scalar
	
	matrix X(4,4);	matrix Y(4,4);	matrix Z(4,4);	matrix W(4,4);	matrix T(4,4);

	matrix  WX(4,4);	matrix  WY(4,4);	matrix  WZ(4,4);	matrix  WT(4,4);	matrix  XY(4,4);
	matrix  XZ(4,4);	matrix  XT(4,4);	matrix  YZ(4,4);	matrix  YT(4,4);	matrix  ZT(4,4);

	matrix  WXY(4,4);	matrix  WXZ(4,4);	matrix  WXT(4,4);	matrix  WYZ(4,4);	matrix  WYT(4,4);
	matrix  WZT(4,4);	matrix  XYZ(4,4);	matrix  XYT(4,4);	matrix  XZT(4,4);	matrix  YZT(4,4);

	matrix  WXYZ(4,4);	matrix  WXYT(4,4);	matrix  WXZT(4,4);	matrix  WYZT(4,4);	matrix  XYZT(4,4);

	matrix  WXYZT(4,4);

/////////////////////////////////////////////////////////////////////

	matrix GA_to_Matrix(GA5_4_1 MV)
{

	int ii,jj;

	matrix spacetime(4,4);
	ex a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A;

	a = MV.q;
	b = MV.w;	c = MV.x;	d = MV.y;	e = MV.z;	f = MV.t;

	g = MV.wx;	h = MV.wy;	j = MV.wz;	k = MV.wt;	l = MV.xy;
	m = MV.xz;	n = MV.xt;	p = MV.yz;	r = MV.yt;	s = MV.zt;

	S = MV.wxy;	R = MV.wxz;	P = MV.wxt;	N = MV.wyz;	M = MV.wyt;
	L = MV.wzt;	K = MV.xyz;	J = MV.xyt;	H = MV.xzt;	G = MV.yzt;

	F = MV.wxyz;	E = MV.wxyt;	D = MV.wxzt;	C = MV.wyzt;	B = MV.xyzt;
	A = MV.wxyzt;


	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			spacetime(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ b*W(ii,jj) + c*X(ii,jj) + d*Y(ii,jj) + e*Z(ii,jj) + f*T(ii,jj)
					+ g*WX(ii,jj) + h*WY(ii,jj) + j*WZ(ii,jj) + k*WT(ii,jj) + l*XY(ii,jj)
					+ m*XZ(ii,jj) + n*XT(ii,jj) + p*YZ(ii,jj) + r*YT(ii,jj) + s*ZT(ii,jj)

					+ S*WXY(ii,jj) + R*WXZ(ii,jj) + P*WXT(ii,jj) + N*WYZ(ii,jj) + M*WYT(ii,jj)
					+ L*WZT(ii,jj) + K*XYZ(ii,jj) + J*XYT(ii,jj) + H*XZT(ii,jj) + G*YZT(ii,jj)
					+ F*WXYZ(ii,jj) + E*WXYT(ii,jj) + D*WXZT(ii,jj) + C*WYZT(ii,jj) + B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}

	return spacetime;
}

/////////////////////////////////////////////////////////////////////

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

	symbol aa("aa"), AA("AA");	// used in complex number format for determinant, etc

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


	int ii,jj;	// loop counters

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

	GA5_4_1 MV1, MV2, MV3, MV4, MV5;
	GA5_4_1 MV6, MV7, MV8, MV9, MV10;

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	cout << "MV1       = " << MV1 << "\n";
/*
	MV5 = Dual(MV1);
	cout << "Dual(MV1) = " << MV5 << "\n\n";

	MV6 = Conjugate(MV1);
	cout << "Conjugate(MV1) = " << MV6 << "\n\n";

	MV7 = Parity(MV1);
	cout << "Parity(MV1) = " << MV7 << "\n\n";

	MV8 = Transpose(MV1);
	cout << "Transpose(MV1) = " << MV8 << "\n\n";

	MV9 = Hermitian(MV1);
	cout << "Hermitian(MV1) = " << MV9 << "\n\n";
	MV9 = MV1*Hermitian(MV1);
	cout << "MV1*Hermitian(MV1) = " << MV9 << "\n\n";

	MV10 = Reverse(MV1);
	cout << "Reverse(MV1) = " << MV10 << "\n\n";

	MV5 = MV1*Reverse(MV1);
	cout << "MV1*Reverse(MV1) = " << MV5 << "\n\n";
*/

////////////////////////////////////  begin matrix calculations


	zero  =  0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0,   0,  0,  0,  0;

	
	Q =  	 1,  0,  0,  0,
		 0,  1,  0,  0,
		 0,  0,  1,  0,
		 0,  0,  0,  1;

	W = 	 0,  0,  0, -I,
		 0,  0,  I,  0,
		 0, -I,  0,  0,
		 I,  0,  0,  0;

	X = 	 1,  0,  0,  0,
		 0, -1,  0,  0,
		 0,  0,  1,  0,
		 0,  0,  0, -1;

	Y = 	 0,  1,  0,  0, 
		 1,  0,  0,  0,  
		 0,  0,  0,  1,
		 0,  0,  1,  0;

	Z = 	 0,  0,  0, -1,
		 0,  0,  1,  0,
		 0,  1,  0,  0, 
		-1,  0,  0,  0;

	T = 	 0, -1,  0,  0,  
		 1,  0,  0,  0,
		 0,  0,  0,  1,
		 0,  0, -1,  0;

//	create our remaining terms

// bivectors
	WX = W.mul(X);	WY = W.mul(Y);	WZ = W.mul(Z);	WT = W.mul(T);	XY = X.mul(Y);	
	XZ = X.mul(Z);	XT = X.mul(T);	YZ = Y.mul(Z);	YT = Y.mul(T);	ZT = Z.mul(T);
 // Trivectors
	WXY = WX.mul(Y);	WXZ = WX.mul(Z);	WXT = WX.mul(T);	WYZ = WY.mul(Z);	WYT = WY.mul(T);
	WZT = WZ.mul(T);	XYZ = XY.mul(Z);	XYT = XY.mul(T);	XZT = XZ.mul(T);	YZT = YZ.mul(T);
// Tetra vectors
	WXYZ = WXY.mul(Z); 	WXYT = WXY.mul(T);	WXZT = WXZ.mul(T);	WYZT = WYZ.mul(T);	XYZT = XYZ.mul(T);
// penta vector
	WXYZT = WXYZ.mul(T);	

//	build generic MV in matrix format

	matrix spacetime(4,4);

	for (ii=0;ii<4;ii++) {
		for (jj=0;jj<4;jj++) {
			spacetime(ii,jj) = expand(
					+ a*Q(ii,jj)
					+ b*W(ii,jj) + c*X(ii,jj) + d*Y(ii,jj) + e*Z(ii,jj) + f*T(ii,jj)
					+ g*WX(ii,jj) + h*WY(ii,jj) + j*WZ(ii,jj) + k*WT(ii,jj) + l*XY(ii,jj)
					+ m*XZ(ii,jj) + n*XT(ii,jj) + p*YZ(ii,jj) + r*YT(ii,jj) + s*ZT(ii,jj)

					+ S*WXY(ii,jj) + R*WXZ(ii,jj) + P*WXT(ii,jj) + N*WYZ(ii,jj) + M*WYT(ii,jj)
					+ L*WZT(ii,jj) + K*XYZ(ii,jj) + J*XYT(ii,jj) + H*XZT(ii,jj) + G*YZT(ii,jj)
					+ F*WXYZ(ii,jj) + E*WXYT(ii,jj) + D*WXZT(ii,jj) + C*WYZT(ii,jj) + B*XYZT(ii,jj)
					+ A*WXYZT(ii,jj) 
                        );
		}
	}

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	ex det1, det2, delta;

//	det1 = spacetime.determinant();
//	det2 = Determinant(MV1);
//	delta = expand(det1 - det2);

//	cout << "spacetime = " << spacetime << "\n";
//	cout << "delta = " << delta << "\n";

	matrix hmm(4,4);

	hmm = GA_to_Matrix(MV1);
	cout << "spacetime = " << spacetime << "\n\n";
	cout << "hmm = " << hmm << "\n\n";
	if(spacetime == hmm) cout << "Matrices agree\n";  else cout << "Clue! - check comparison methods!\n";

//	now we check the transpose function

	matrix argh(4,4);
	argh = spacetime.transpose();
	hmm = GA_to_Matrix(Transpose(MV1));
	if (hmm == argh) cout << "Transpose Matrices agree\n";  else cout << "Clue! - Transpose failure!\n";

// check assignment operator problems

	argh = spacetime.transpose();
	MV2 = Transpose(MV1);
	hmm = GA_to_Matrix(MV2);
	if (hmm == argh) cout << "Transpose Matrices agree\n";  else cout << "Clue! - Transpose failure!\n";

// check determinant again

	det1 = Determinant(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	cout << "delta = " << delta << "\n";
	if(det1 == det2) cout << "Determinants agree \n"; else cout << "must expand before comparisons\n";

// now I want to look at the complex determinant for a prime number test vector for the failed unary operations

	MV2 = OverBar(MV1) ;  
	cout << "OverBar \n" << MV2 << "\n\n";

	MV2 = UnderBar(MV1) ;  
	cout << "UnderBar \n" << MV2 << "\n\n";

	MV2 = Reverse(MV1) ;  
	cout << "Reverse \n" << MV2 << "\n\n";

//	MV2 = Reverse_Vector_Quad(MV1) ;  
//	cout << "Reverse_Vector_Quad \n" << MV2 << "\n\n";

	MV2 = Involution(MV1) ;  
	cout << "Involution \n" << MV2 << "\n\n";

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

	MV2 = Involution(MV1) ;  
	cout << "Involution \n" << Determinant(MV2) << "\n\n";

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
	 

	MV3 = MV1*OverBar(MV1) ;
	cout << "MV1*OverBar(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*UnderBar(MV1) ;
	cout << "MV1*UnderBar(MV1) = " << MV3 << "\n";
	cout << Determinant(MV3) << "\n\n";

	MV3 = MV1*Involution(MV1) ;
	cout << "MV1*Involution(MV1) = " << MV3 << "\n";
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

	GA5_4_1 Blade[32], Zero;
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

	return 0;
}

/*
GA5_4_1 Parity(const GA5_4_1 &a)	// Corresponds to parity transform: w -> -w, x -> -x, y -> -y, z ->-z, t -> -t
GA5_4_1 Dual(const GA5_4_1 &a)	// multiply by pseudoscalar. LeftDual = RightDual = Dual
GA5_4_1 ComplexConjugate(const GA5_4_1 &a)	// negate every term including w
GA5_4_1 Transpose(const GA5_4_1 &a)	
GA5_4_1 Hermitian(const GA5_4_1 &a)	
GA5_4_1 Reverse(const GA5_4_1 &a)
	r = Set(  3,    5,  7, 11, 13, 17,
		 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 
		 19, 23, 29, 31, 37,    41) ; 
	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

*/






