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

int main(void)
{

	symbol a_q("a");
	symbol a_x("b") ,  a_y("c") ,  a_z("d"),   a_t("e");
	symbol a_xy("f"),  a_xz("g"),  a_yz("h"),  a_xt("i"), a_yt("j"), a_zt("k");
	symbol a_xyz("l"), a_xyt("m"), a_xzt("n"), a_yzt("o");
	symbol a_xyzt("p");

	Mink r,s,t,u,w;
	Mink Zero;   Zero = Mink();

	ex Det1, Det2;

	r = Mink(a_q,a_x,a_y,a_z,a_t,a_xy,a_xz,a_yz,a_xt,a_yt,a_zt,a_xyz,a_xyt,a_xzt,a_yzt,a_xyzt);  
	s  = Reverse(r);

	w = Product(r,s) - Product(s,r);

	if(w == Zero) cout << "r and Reverse(r) commute\n"; else cout << "r and Reverse(r) do not commute\n";
	cout << "Commutator = " << w << "\n";

	Det1 = Determinant(r);

	Mink A, B, C;
	ex a, b, c, d, e, f, delta;

	A = r;

	B = Reverse(A);
	C = B*A;

	a = C.q;
	b = C.x;
	c = C.y;
	d = C.z;
	e = C.t;

	f = C.xyzt;
	Det2 = expand(a*a - b*b - c*c - d*d + e*e + f*f); // works due to square of component terms being used

	if (Det1 == Det2) cout << "Determinants agree\n"; else cout << "Determinants disagree\n";
	delta = expand(Det1 - Det2);
	cout << "Delta = " << delta << "\n";

	cout << "FYI I*I = " << I*I << "\n";

	return 0;
}

