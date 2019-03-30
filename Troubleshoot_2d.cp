// compare matrix product and geometric product

// compile by
// g++ Demo_GA2E.cp -l ginac -l cln


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>

#include "GA2E_Routines.cp"

using namespace std;
using namespace GiNaC;


///////////////////////////////////

int main(void)
{

	
//	int i,j,k;

	GA2E check;
	GA2E Zero;   Zero = GA2E();

	GA2E q, x,y, xy;
	GA2E r,s,t,u;

	symbol a_q("a"), a_x("b"), a_y("c"), a_xy("d") ;
	symbol b_q("A"), b_x("B") ,  b_y("C"), b_xy("D") ;
	symbol c_q("c.a"), c_x("c.b"), c_y("c.c"), c_xy("c.d") ;
	symbol d_q("d.a"), d_x("d.b"), d_y("d.c"), d_xy("d.d");

////////////////////////////////////

	q  = GA2E(1, 0,0, 0);
	x  = GA2E(0, 1,0, 0);
	y  = GA2E(0, 0,1, 0);
	xy = GA2E(0, 0,0, 1);

// Check the OverBar and UnderBar functions

// OverBar // a blade wedge b blade = pseudovector blade xy

	r = Wedge( q, OverBar( q));  cout << "r = " << r << " expect (0, 0,0, 1) \n";
	r = Wedge( x, OverBar( x));  cout << "r = " << r << " expect (0, 0,0, 1) \n";
	r = Wedge( y, OverBar( y));  cout << "r = " << r << " expect (0, 0,0, 1) \n";
	r = Wedge(xy, OverBar(xy));  cout << "r = " << r << " expect (0, 0,0, 1) \n\n";
// UnderBar// b blade wedge a blade = pseudovector blade
	r = Wedge(UnderBar( q),  q);  cout << "r = " << r << " expect (0, 0,0, 1) \n";
	r = Wedge(UnderBar( x),  x);  cout << "r = " << r << " expect (0, 0,0, 1) \n";
	r = Wedge(UnderBar( y),  y);  cout << "r = " << r << " expect (0, 0,0, 1) \n";
	r = Wedge(UnderBar(xy), xy);  cout << "r = " << r << " expect (0, 0,0, 1) \n\n";

/////////////////////////////////  

	r = GA2E( 3,  5,  7, 11) ; 
	s = OverBar(r);
	t = Product(r,s);
	cout << "r = " << r << "\n";
	cout << "s = " << s << "\n";
	cout << "t = r*OverBar(r) = " << t << "\n\n";

	r = GA2E( 3,  5,  7, 11) ; 
	s = UnderBar(r);
	t = Product(r,s);
	cout << "r = " << r << "\n";
	cout << "s = " << s << "\n";
	cout << "t = r*UnderBar(r) = " << t << "\n\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = OverBar(r);
	t = Product(r,s);
	cout << "r = " << r << "\n";
	cout << "s = " << s << "\n";
	cout << "t = r*OverBar(r) = " << t << "\n\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = UnderBar(r);
	t = Product(r,s);
	cout << "r = " << r << "\n";
	cout << "s = " << s << "\n";
	cout << "t = r*UnderBar(r) = " << t << "\n\n";

	cout << "Determinant(r) = " << Determinant(r) << "\n";
	cout << "Determinant(r) - t.xy*t.xy = " << expand( Determinant(r) - t.xy*t.xy) << "\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = OverBar(r);
	if(Determinant(r) == Determinant(s)) cout << "OverBar conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = UnderBar(r);
	if(Determinant(r) == Determinant(s)) cout << "UnderBar conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = Dual(r);
	if(Determinant(r) == Determinant(s)) cout << "Dual conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = DorstDual(r);
	if(Determinant(r) == Determinant(s)) cout << "DorstDual conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = DorstUnDual(r);
	if(Determinant(r) == Determinant(s)) cout << "DorstUnDual conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = Conjugation(r);
	if(Determinant(r) == Determinant(s)) cout << "Conjugation conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = CliffordConjugation(r);
	if(Determinant(r) == Determinant(s)) cout << "CliffordConjugation conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = Transpose(r);
	if(Determinant(r) == Determinant(s)) cout << "Transpose conserves the determinant\n";

	r = GA2E(a_q, a_x, a_y, a_xy) ; 
	s = Reciprocal(r);
	t = Product(r,s);
	cout << "1/r = " << s << "\n";
	cout << "(r*(1/r)) = " << t << " expect (1, 0,0, 0); \n";

/////////////////////////////////  

	return 0;
}



