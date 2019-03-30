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

/////////////////////////////////// Test Determinant ////////////////////////////////


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

	ex det1, det2, delta;
/*
	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	det1 = Determinant(MV1);

	MV2 = Transpose(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Transpose delta = " << delta << "  (expect 0)\n\n";  // good
	else cout << "Transpose does not conserve determinant\n";

	MV2 = ComplexConjugate(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "ComplexConjugate delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "ComplexConjugate does not conserve determinant\n";

	MV2 = Parity(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Parity delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "Parity does not conserve determinant\n";

	MV2 = Reverse(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Reverse delta = " << delta << "  (expect 0)\n\n";  // works
	else cout << "Reverse does not conserve determinant\n";

	MV2 = CliffordConjugation(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "CliffordConjugation delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "CliffordConjugation does not conserve determinant\n";

	MV2 = Hermitian(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Hermitian delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "Hermitian does not conserve determinant\n";
*/
	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	det1 = Determinant(MV1);

	MV2 = Transpose(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Transpose delta = " << delta << "  (expect 0)\n\n";  // good
	else cout << "Transpose does not conserve determinant\ndet1 = " << det1 << "\ndet2 = " << det2 << "\n\n";

	MV2 = ComplexConjugate(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "ComplexConjugate delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "ComplexConjugate does not conserve determinant\ndet1 = " << det1 << "\ndet2 = " << det2 << "\n\n";

	MV2 = Parity(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Parity delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "Parity does not conserve determinant\ndet1 = " << det1 << "\ndet2 = " << det2 << "\n\n";

	MV2 = Reverse(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Reverse delta = " << delta << "  (expect 0)\n\n";  // works
	else cout << "Reverse does not conserve determinant\ndet1 = " << det1 << "\ndet2 = " << det2 << "\n\n";

	MV2 = CliffordConjugation(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "CliffordConjugation delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "CliffordConjugation does not conserve determinant\ndet1 = " << det1 << "\ndet2 = " << det2 << "\n\n";

	MV2 = Hermitian(MV1);
	det2 = Determinant(MV2);
	delta = expand(det1 - det2);
	if (delta == 0) cout << "Hermitian delta = " << delta << "  (expect 0)\n\n";  // fails
	else cout << "Hermitian does not conserve determinant\ndet1 = " << det1 << "\ndet2 = " << det2 << "\n\n";


////////////////////////////////////
/*
//	Verify library determinant versus matrix


	matrix Zero(4,4);
	matrix q(4,4), x(4,4),y(4,4),z(4,4), xy(4,4),yz(4,4),xz(4,4), xyz(4,4);
	matrix t(4,4), xt(4,4),yt(4,4),zt(4,4), xyt(4,4),yzt(4,4),xzt(4,4), xyzt(4,4);
	matrix w(4,4), wx(4,4),wy(4,4),wz(4,4), wxy(4,4),wyz(4,4),wxz(4,4), wxyz(4,4);
	matrix wt(4,4), wxt(4,4),wyt(4,4),wzt(4,4), wxyt(4,4),wyzt(4,4),wxzt(4,4), wxyzt(4,4);

	matrix MV(4,4);


	Zero = 	 0, 0, 0, 0,
		 0, 0, 0, 0,
		 0, 0, 0, 0,
		 0, 0, 0, 0  ;

	q = 	 1, 0, 0, 0,
		 0, 1, 0, 0,
		 0, 0, 1, 0,
		 0, 0, 0, 1  ;

	x = 	 0, 0, 1, 0,
		 0, 0, 0, 1,
		 1, 0, 0, 0,
		 0, 1, 0, 0  ;

	y = 	 1, 0, 0, 0,
		 0, 1, 0, 0,
		 0, 0,-1, 0,
		 0, 0, 0,-1  ;

	z = 	 0, 0, 0, 1,
		 0, 0,-1, 0,
		 0,-1, 0, 0,
		 1, 0, 0, 0  ;

	t = 	 0, 0,-1, 0,
		 0, 0, 0, 1,
		 1, 0, 0, 0,
		 0,-1, 0, 0  ;

	wxyzt =  I, 0, 0, 0,
		 0, I, 0, 0,
		 0, 0, I, 0,
		 0, 0, 0, I  ;

//////////////////////////////////////////////////////

//   calculate remaining terms. Synthesize up for xyzt

	xy = x.mul(y);
	xz = x.mul(z);
	yz = y.mul(z);
	xyz = xy.mul(z);

	xt = x.mul(t);
	yt = y.mul(t);
	zt = z.mul(t);

	xyt = xy.mul(t);
	xzt = xz.mul(t);
	yzt = yz.mul(t);
	xyzt = xyz.mul(t);

//  I assigned wxyzt to be I. Synthesize down for the remaining w terms

	wxyz = Zero.sub(wxyzt.mul(t));	// t^2 = -1
	wxyt = Zero.sub(wxyzt.mul(z));
	wxzt = wxyzt.mul(y);
	wyzt = Zero.sub(wxyzt.mul(x));

	wxy = wxyz.mul(z);
	wxz = Zero.sub(wxyz.mul(y));
	wxt = Zero.sub(wxyt.mul(y));
	wyz = Zero.sub(wyzt.mul(t));
	wyt = wxyt.mul(x);
	wzt = wyzt.mul(y);

	wx = wxy.mul(y);
	wy = Zero.sub(wxy.mul(x));
	wz = Zero.sub(wxz.mul(x));
	wt = Zero.sub(wxt.mul(x));

	w = wx.mul(x);

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	det1 = Determinant(MV1);

//	ex q,   w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt,
//		wxy,wxz,wxt,wyz,wyt,wzt,xyz,xyt,xzt,yzt,  
//		wxyz,wxyt,wxzt,wyzt,xyzt,  wxyzt;

	MV = Zero;
	MV = MV.add(q.mul_scalar(a));

	MV = MV.add(w.mul_scalar(b));
	MV = MV.add(x.mul_scalar(c));
	MV = MV.add(y.mul_scalar(d));
	MV = MV.add(z.mul_scalar(e));
	MV = MV.add(t.mul_scalar(f));

	MV = MV.add(wx.mul_scalar(g));
	MV = MV.add(wy.mul_scalar(h));
	MV = MV.add(wz.mul_scalar(j));
	MV = MV.add(wt.mul_scalar(k));
	MV = MV.add(xy.mul_scalar(l));
	MV = MV.add(xz.mul_scalar(m));
	MV = MV.add(xt.mul_scalar(n));
	MV = MV.add(yz.mul_scalar(p));
	MV = MV.add(yt.mul_scalar(r));
	MV = MV.add(zt.mul_scalar(s));

	MV = MV.add(wxy.mul_scalar(S));
	MV = MV.add(wxz.mul_scalar(R));
	MV = MV.add(wxt.mul_scalar(P));
	MV = MV.add(wyz.mul_scalar(N));
	MV = MV.add(wyt.mul_scalar(M));
	MV = MV.add(wzt.mul_scalar(L));
	MV = MV.add(xyz.mul_scalar(K));
	MV = MV.add(xyt.mul_scalar(J));
	MV = MV.add(xzt.mul_scalar(H));
	MV = MV.add(yzt.mul_scalar(G));

	MV = MV.add(wxyz.mul_scalar(F));
	MV = MV.add(wxyt.mul_scalar(E));
	MV = MV.add(wxzt.mul_scalar(D));
	MV = MV.add(wyzt.mul_scalar(C));
	MV = MV.add(xyzt.mul_scalar(B));

	MV = MV.add(wxyzt.mul_scalar(A));

	det2 = MV.determinant();
	delta = expand(det1 - det2);
	cout << "\n\nComparing matrix versus GA determinant formulas . . .\n";
	if (delta == 0) cout << "Determinants Agree\n\n";  // Determinants agree. Good.
	else cout << "Determinants differ. \ndet1 = " << det1 << "\ndet2 = " << det2 << "\n\n";

////////////////////////////////////

//	See if the differences are always conjugation
	ex sum;

	cout << "\nCut line here ----------------------------------\n\n";

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);
	det1 = Determinant(MV1);

	MV2 = ComplexConjugate(MV1);
	det2 = Determinant(MV2);	// conjugate indeed
	sum = expand((det1 + det2)/2);  // see if imaginary component is zero. Expect no letter I in text result
	cout << "\n\nComplexConjugate sum is " << sum << "\n\n";

	MV2 = Parity(MV1);
	det2 = Determinant(MV2);	// conjugate indeed
	sum = expand((det1 + det2)/2);  // see if imaginary component is zero
	cout << "\n\nParity sum is " << sum << "\n\n";

	MV2 = CliffordConjugation(MV1);
	det2 = Determinant(MV2);	// conjugate indeed
	sum = expand((det1 + det2)/2);  // see if imaginary component is zero
	cout << "\n\nCliffordConjugation sum is " << sum << "\n\n";

	MV2 = Hermitian(MV1);
	det2 = Determinant(MV2);	// conjugate indeed
	sum = expand((det1 + det2)/2);  // see if imaginary component is zero
	cout << "\n\nHermitian sum is " << sum << "\n\n";
*/
////////////////////////////////////

	ex Sum_of_Squares_1, Sum_of_Squares_2;

	MV1 = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 
	det1 = Determinant(MV1);
	Sum_of_Squares_1 = real_part(det1)*real_part(det1) + imag_part(det1)*imag_part(det1);

	MV2 = Reverse(MV1);
	det2 = Determinant(MV2);
	cout << "Reverse : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

	MV2 = Reverse_Vector_Quad(MV1);
	det2 = Determinant(MV2);
	cout << "Reverse_Vector_Quad : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

	MV2 = Conjugation(MV1);
	det2 = Determinant(MV2);
	cout << "Conjugation : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

	MV2 = CliffordConjugation(MV1);
	det2 = Determinant(MV2);
	cout << "CliffordConjugation : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

	MV2 = Parity(MV1);
	det2 = Determinant(MV2);
	cout << "Parity : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

	MV2 = ComplexConjugate(MV1);
	det2 = Determinant(MV2);
	cout << "ComplexConjugate : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

	MV2 = Hermitian(MV1);
	det2 = Determinant(MV2);
	cout << "Hermitian : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

	MV2 = Conjugate(MV1);
	det2 = Determinant(MV2);
	cout << "Conjugate : " << det1 << "    " << det2 << "\n";
	Sum_of_Squares_2 = real_part(det2)*real_part(det2) + imag_part(det2)*imag_part(det2);
	if(Sum_of_Squares_1 == Sum_of_Squares_2) cout << "Magnitudes agree\n";

////////////////////////////////////



	return 0;
}

