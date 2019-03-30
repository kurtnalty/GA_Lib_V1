#include <stdio.h>
#include <stdlib.h>

#include "Mink_Routines.c"

/////////////////////////////////// Demonstrate Mink Routines ////////////////////////////////


int main(void)
{
	
	Mink r,s,t,u;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Mink_Set(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

	double det_ref, det;
	det_ref = Mink_Determinant(r);
	printf("\nMink Routines Demo\n\n");

//////////////////////////////////////////////////////

// void Mink_PrintlnMV(Mink v) 

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("Mink_PrintlnMV Demo\nr = ");  Mink_PrintlnMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void Mink_PrintMV(Mink v) 

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("Mink_PrintMV Demo\nr = ");  Mink_PrintMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_OverBar(Mink a) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_OverBar(r);

	printf("Mink_OverBar Demo\n");
	printf("r =              ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_OverBar(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =       ");  Mink_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_UnderBar(Mink a) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_UnderBar(r);

	printf("Mink_UnderBar Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_UnderBar(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	t = Mink_OverBar(Mink_UnderBar(r));
	printf("t = Mink_OverBar(Mink_UnderBar(r)) =        ");  Mink_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Reverse(Mink w) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Reverse(r);

	printf("Mink_Reverse Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Reverse(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	t = Mink_Reverse(Mink_Reverse(r));
	printf("t = Mink_Reverse(Mink_Reverse(r)) =        ");  Mink_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Parity(Mink w) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Parity(r);

	printf("Mink_Parity Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Parity(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	t = Mink_Parity(Mink_Parity(r));
	printf("t = Mink_Parity(Mink_Parity(r)) =        ");  Mink_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Transpose(Mink w) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Transpose(r);

	printf("Mink_Transpose Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Transpose(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	t = Mink_Transpose(Mink_Transpose(r));
	printf("t = Mink_Transpose(Mink_Transpose(r)) =        ");  Mink_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Conjugation(Mink w) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Conjugation(r);

	printf("Mink_Conjugation Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Conjugation(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	t = Mink_Conjugation(Mink_Conjugation(r));
	printf("t = Mink_Conjugation(Mink_Conjugation(r)) =        ");  Mink_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Clifford_Conjugation(Mink w) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Clifford_Conjugation(r);

	printf("Mink_Clifford_Conjugation Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Clifford_Conjugation(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	t = Mink_Clifford_Conjugation(Mink_Clifford_Conjugation(r));
	printf("t = Mink_Clifford_Conjugation(Mink_Clifford_Conjugation(r)) =        ");  Mink_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Dual(Mink w) ; 

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Dual(r);

	printf("Mink_Dual Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Dual(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	t = Mink_Dual(Mink_Dual(r));
	printf("t = Mink_Dual(Mink_Dual(r)) =        ");  Mink_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Dorst_Dual(Mink a) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Dorst_Dual(r);

	printf("Mink_Dorst_Dual Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Reverse(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Dorst_UnDual(Mink a) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = Mink_Dorst_UnDual(r);

	printf("Mink_Dorst_UnDual Demo\n");
	printf("r =               ");  Mink_PrintMV(r);  printf("\n");
	printf("u = Mink_Dorst_UnDual(r) = ");  Mink_PrintMV(u);  printf("\n");
	det = Mink_Determinant(u);
	if(det == det_ref) printf("Mink_Determinant is conserved\n"); else printf("Mink_Determinant is not conserved\n");
	s = Mink_Product(r,u);
	printf("s = r*u =        ");  Mink_PrintMV(s);  printf("\n");

	s = Mink_Dorst_Dual(r);
	u = Mink_Dorst_UnDual(s);

	printf("Mink_Dorst_UnDual(Mink_Dorst_Dual(r)) = "); Mink_PrintMV(u); printf("\n");

	s = Mink_Dorst_UnDual(r);
	u = Mink_Dorst_Dual(s);

	printf("Mink_Dorst_Dual(Mink_Dorst_UnDual(r)) = "); Mink_PrintMV(u); printf("\n\n");


//////////////////////////////////////////////////////

// Mink Mink_Add(Mink u,Mink v)

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	u = Mink_Add(r, s);

	printf("  r = ");  Mink_PrintMV(r);  printf("\n");
	printf("  s = ");  Mink_PrintMV(s);  printf("\n");
	printf("r+s = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Subtract(Mink u, Mink v)

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	u = Mink_Subtract(r, s);

	printf("  r = ");  Mink_PrintMV(r);  printf("\n");
	printf("  s = ");  Mink_PrintMV(s);  printf("\n");
	printf("r-s = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int Mink_Equal(Mink u, Mink v) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 

	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	printf("t = ");  Mink_PrintMV(t);  printf("\n");
	if (Mink_Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (Mink_Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int Mink_Not_Equal(Mink u, Mink v) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 

	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	printf("t = ");  Mink_PrintMV(t);  printf("\n");
	if (Mink_Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (Mink_Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// Mink Mink_Product(Mink a, Mink b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Product(r,s);
	printf("u = r*s = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Divide_By_Constant(Mink u, double a) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	u = Mink_Divide_By_Constant(r,2);
	printf("u = r/2 = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Wedge(Mink a, Mink b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Wedge(r,s);
	printf("u = Mink_Wedge(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_AntiWedge(Mink a, Mink b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_AntiWedge(r,s);
	printf("u = Mink_AntiWedge(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Regressive(Mink a, Mink b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Regressive(r,s);
	printf("u = Mink_Regressive(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Regressive_Via_Formula(Mink a, Mink b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Regressive_Via_Formula(r,s);
	printf("u = Mink_Regressive_Via_Formula(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Lower_Right_Via_Formula(Mink a, Mink b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Lower_Right_Via_Formula(r,s);
	printf("u = Mink_Lower_Right_Via_Formula(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Expander(const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Expander(r,s);
	printf("u = Mink_Expander(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Conserver(const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Conserver(r,s);
	printf("u = Mink_Conserver(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Shrinker(const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Shrinker(r,s);
	printf("u = Mink_Shrinker(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// Mink Mink_Symmetric(const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Symmetric(r,s);
	printf("u = Mink_Symmetric(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_AntiSymmetric(const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_AntiSymmetric(r,s);
	printf("u = Mink_AntiSymmetric(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Inner(const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Inner(r,s);
	printf("u = Mink_Inner(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Left_Contraction (const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Left_Contraction(r,s);
	printf("u = Mink_Left_Contraction(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Right_Contraction (const Mink &a, const Mink &b) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	printf("s = ");  Mink_PrintMV(s);  printf("\n");
	u = Mink_Right_Contraction(r,s);
	printf("u = Mink_Right_Contraction(r, s) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex Mink_Determinant(Mink A) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	det = Mink_Determinant(r);
	printf("u = Mink_Determinant(r) = %10.3e \n\n",det); 


//////////////////////////////////////////////////////

// Mink Mink_Adjugate(Mink V) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	u = Mink_Adjugate(r);
	printf("u = Mink_Adjugate(r) = ");  Mink_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// Mink Mink_Reciprocal(Mink a) ;

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  Mink_PrintMV(r);  printf("\n");
	u = Mink_Reciprocal(r);
	printf("u = Mink_Reciprocal(r) = ");  Mink_PrintMV(u);  printf("\n\n");
	s = Mink_Product(r,u);
	printf("s = r*(1/r) = ");  Mink_PrintMV(s);  printf("\n\n");

/////////////////////////////////  

// Verify 64 Comps preserving the determinant

	int i;
	Mink X, Y;

	X = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	det_ref = Mink_Determinant(X);

	printf("Verifying 64 Mink_Determinant preserving Mink_Comp transforms\n");
	for (i=0; i<64; i++) {
		Y = Mink_Comp(X,i);
		det = Mink_Determinant(Y);
		if (det == det_ref) {
			printf("+");
		}
		else printf("Mink_Comp determinant mismatch at index %d \n", i);
	}
	printf("\n");


// Verify 64 Comps preserving the determinant

	X = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	det_ref = Mink_Determinant(X);

	printf("Verifying 72 Mink_Determinant preserving Mink_Magic transforms\n");
	for (i=0; i<72; i++) {
		Y = Mink_Magic(X,i);
		det = Mink_Determinant(Y);
		if (det == det_ref) {
			printf("+");
		}
		else printf("Mink_Magic determinant mismatch at index %d \n", i);
	}
	printf("\n");

/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	r = Mink_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = Mink_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = Mink_Set(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

*/





