#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

#include "Quadratic_Cubic_Quartic.cp"

// ************************************************************************

int main(void)
{
	ex r1, r2, r3, r4;
	ex R1, R2, R3, R4;
	ex A,B,C,D,E; 

	R1 = 1 + 2*I;
	R2 = 3 + 4*I;
	R3 = 5 + 6*I;
	R4 = 7 + 8*I;

	Build_Quadratic(&A, &B, &C, R1, R2);
	cout << "A = " << A << "\n";
	cout << "B = " << B << "\n";
	cout << "C = " << C << "\n";

	Solve_Quadratic(A, B, C, &r1, &r2);
	cout << "r1 = " << r1 << "\n";
	cout << "r2 = " << r2 << "\n\n";


	Build_Cubic(&A, &B, &C, &D, R1, R2, R3);
	cout << "A = " << A << "\n";
	cout << "B = " << B << "\n";
	cout << "C = " << C << "\n";
	cout << "D = " << D << "\n";

	Solve_Cubic(A, B, C, D, &r1, &r2, &r3);
	cout << "r1 = " << r1 << "\n";
	cout << "r2 = " << r2 << "\n";
	cout << "r3 = " << r3 << "\n\n";


	Build_Quartic(&A, &B, &C, &D, &E, R1, R2, R3, R4);
	cout << "A = " << A << "\n";
	cout << "B = " << B << "\n";
	cout << "C = " << C << "\n";
	cout << "D = " << D << "\n";
	cout << "E = " << E << "\n";

	Solve_Quartic(A, B, C, D, E, &r1, &r2, &r3, &r4);
	cout << "r1 = " << r1 << "\n";
	cout << "r2 = " << r2 << "\n";
	cout << "r3 = " << r3 << "\n";
	cout << "r4 = " << r4 << "\n\n";


/////////////////////////////////////////

	return 0;
}



