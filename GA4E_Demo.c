#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GA4E_Routines.c"

/////////////////////////////////// Demonstrate GA4E Routines ////////////////////////////////


int main(void)
{
	
	GA4E r,s,t,u;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = GA4E_Set(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

	double det_ref, det;
	det_ref = GA4E_Determinant(r);
	printf("\nGA4E Routines Demo\n\n");

//////////////////////////////////////////////////////

// void GA4E_PrintlnMV(GA4E v) 

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("GA4E_PrintlnMV Demo\nr = ");  GA4E_PrintlnMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void GA4E_PrintMV(GA4E v) 

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("GA4E_PrintMV Demo\nr = ");  GA4E_PrintMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_OverBar(GA4E a) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_OverBar(r);

	printf("GA4E_OverBar Demo\n");
	printf("r =              ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_OverBar(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =       ");  GA4E_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_UnderBar(GA4E a) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_UnderBar(r);

	printf("GA4E_UnderBar Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_UnderBar(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	t = GA4E_OverBar(GA4E_UnderBar(r));
	printf("t = GA4E_OverBar(GA4E_UnderBar(r)) =        ");  GA4E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Reverse(GA4E w) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_Reverse(r);

	printf("GA4E_Reverse Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_Reverse(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	t = GA4E_Reverse(GA4E_Reverse(r));
	printf("t = GA4E_Reverse(GA4E_Reverse(r)) =        ");  GA4E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Involution(GA4E w) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_Involution(r);

	printf("GA4E_Involution Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_Involution(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	t = GA4E_Involution(GA4E_Involution(r));
	printf("t = GA4E_Involution(GA4E_Involution(r)) =        ");  GA4E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Transpose(GA4E w) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_Transpose(r);

	printf("GA4E_Transpose Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_Transpose(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	t = GA4E_Transpose(GA4E_Transpose(r));
	printf("t = GA4E_Transpose(GA4E_Transpose(r)) =        ");  GA4E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Conjugation(GA4E w) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_Conjugation(r);

	printf("GA4E_Conjugation Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_Conjugation(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	t = GA4E_Conjugation(GA4E_Conjugation(r));
	printf("t = GA4E_Conjugation(GA4E_Conjugation(r)) =        ");  GA4E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_CliffordConjugation(GA4E w) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_CliffordConjugation(r);

	printf("GA4E_CliffordConjugation Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_CliffordConjugation(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	t = GA4E_CliffordConjugation(GA4E_CliffordConjugation(r));
	printf("t = GA4E_CliffordConjugation(GA4E_CliffordConjugation(r)) =        ");  GA4E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Dual(GA4E w) ; 

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_Dual(r);

	printf("GA4E_Dual Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_Dual(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	t = GA4E_Dual(GA4E_Dual(r));
	printf("t = GA4E_Dual(GA4E_Dual(r)) =        ");  GA4E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_DorstDual(GA4E a) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_DorstDual(r);

	printf("GA4E_DorstDual Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_Reverse(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_DorstUnDual(GA4E a) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	u = GA4E_DorstUnDual(r);

	printf("GA4E_DorstUnDual Demo\n");
	printf("r =               ");  GA4E_PrintMV(r);  printf("\n");
	printf("u = GA4E_DorstUnDual(r) = ");  GA4E_PrintMV(u);  printf("\n");
	det = GA4E_Determinant(u);
	if(det == det_ref) printf("GA4E_Determinant is conserved\n"); else printf("GA4E_Determinant is not conserved\n");
	s = GA4E_Product(r,u);
	printf("s = r*u =        ");  GA4E_PrintMV(s);  printf("\n");

	s = GA4E_DorstDual(r);
	u = GA4E_DorstUnDual(s);

	printf("GA4E_DorstUnDual(GA4E_DorstDual(r)) = "); GA4E_PrintMV(u); printf("\n");

	s = GA4E_DorstUnDual(r);
	u = GA4E_DorstDual(s);

	printf("GA4E_DorstDual(GA4E_DorstUnDual(r)) = "); GA4E_PrintMV(u); printf("\n\n");


//////////////////////////////////////////////////////

// GA4E GA4E_Add(GA4E u,GA4E v)

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	u = GA4E_Add(r, s);

	printf("  r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("  s = ");  GA4E_PrintMV(s);  printf("\n");
	printf("r+s = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Subtract(GA4E u, GA4E v)

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	u = GA4E_Subtract(r, s);

	printf("  r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("  s = ");  GA4E_PrintMV(s);  printf("\n");
	printf("r-s = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int GA4E_Equal(GA4E u, GA4E v) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 

	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	printf("t = ");  GA4E_PrintMV(t);  printf("\n");
	if (GA4E_Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (GA4E_Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int GA4E_Not_Equal(GA4E u, GA4E v) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 

	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	printf("t = ");  GA4E_PrintMV(t);  printf("\n");
	if (GA4E_Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (GA4E_Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Product(GA4E a, GA4E b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Product(r,s);
	printf("u = r*s = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Divide_By_Constant(GA4E u, double a) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	u = GA4E_Divide_By_Constant(r,2);
	printf("u = r/2 = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Wedge(GA4E a, GA4E b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Wedge(r,s);
	printf("u = GA4E_Wedge(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_AntiWedge(GA4E a, GA4E b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_AntiWedge(r,s);
	printf("u = GA4E_AntiWedge(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Regressive(GA4E a, GA4E b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Regressive(r,s);
	printf("u = GA4E_Regressive(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_RegressiveViaFormula(GA4E a, GA4E b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_RegressiveViaFormula(r,s);
	printf("u = GA4E_RegressiveViaFormula(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_LowerRightViaFormula(GA4E a, GA4E b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_LowerRightViaFormula(r,s);
	printf("u = GA4E_LowerRightViaFormula(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Expander(const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Expander(r,s);
	printf("u = GA4E_Expander(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Conserver(const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Conserver(r,s);
	printf("u = GA4E_Conserver(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Shrinker(const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Shrinker(r,s);
	printf("u = GA4E_Shrinker(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// GA4E GA4E_Symmetric(const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Symmetric(r,s);
	printf("u = GA4E_Symmetric(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_AntiSymmetric(const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_AntiSymmetric(r,s);
	printf("u = GA4E_AntiSymmetric(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Inner(const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_Inner(r,s);
	printf("u = GA4E_Inner(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_LeftContraction (const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_LeftContraction(r,s);
	printf("u = GA4E_LeftContraction(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_RightContraction (const GA4E &a, const GA4E &b) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	printf("s = ");  GA4E_PrintMV(s);  printf("\n");
	u = GA4E_RightContraction(r,s);
	printf("u = GA4E_RightContraction(r, s) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex GA4E_Determinant(GA4E A) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	det = GA4E_Determinant(r);
	printf("u = GA4E_Determinant(r) = %10.3e \n\n",det); 


//////////////////////////////////////////////////////

// GA4E GA4E_Adjugate(GA4E V) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	u = GA4E_Adjugate(r);
	printf("u = GA4E_Adjugate(r) = ");  GA4E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA4E GA4E_Reciprocal(GA4E a) ;

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	printf("r = ");  GA4E_PrintMV(r);  printf("\n");
	u = GA4E_Reciprocal(r);
	printf("u = GA4E_Reciprocal(r) = ");  GA4E_PrintMV(u);  printf("\n\n");
	s = GA4E_Product(r,u);
	printf("s = r*(1/r) = ");  GA4E_PrintMV(s);  printf("\n\n");

/////////////////////////////////  

// Verify 64 GA4E_Comps preserving the determinant

	int i;
	GA4E X, Y;

	X = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	det_ref = GA4E_Determinant(X);

	printf("Verifying 64 GA4E_Determinant preserving Comp transforms\n");
	for (i=0; i<64; i++) {
		Y = GA4E_Comp(X,i);
		det = GA4E_Determinant(Y);
		if (det == det_ref) {
			printf("+");
		}
		else printf("Comp determinant mismatch at index %d \n", i);
	}
	printf("\n");


// Verify 120 GA4E_Magic preserving the determinant

	X = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	det_ref = GA4E_Determinant(X);

	printf("Verifying 120 GA4E_Determinant preserving GA4E_Magic transforms\n");
	for (i=0; i<120; i++) {
		Y = GA4E_Magic(X,i);
		det = GA4E_Determinant(Y);
		if (det == det_ref) {
			printf("+");
		}
		else printf("GA4E_Magic determinant mismatch at index %d \n", i);
	}
	printf("\n");


/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	r = GA4E_Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 
	s = GA4E_Set( 61,   67, 71, 73, 79,     83, 89, 97,101,103,107,   109,113,127,131,   137) ; 
	t = GA4E_Set(139,  149,151,157,163,    167,173,179,181,191,193,   197,199,211,223,   227) ; 

*/





