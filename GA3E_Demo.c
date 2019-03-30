#include <stdio.h>
#include <stdlib.h>

#include "GA3E_Routines.c"

/////////////////////////////////// Demonstrate GA3E Routines ////////////////////////////////


int main(void)
{
	

	GA3E r,s,t,u;

//	GA3E X, Y, Z;

//	X = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
//	Y = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
//	Z = GA3E_Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
//	t = GA3E_Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

	double det_ref, det;
	det_ref = GA3E_Determinant(r);
	printf("\nGA3E Routines Demo\n\n");

//////////////////////////////////////////////////////

// void GA3E_PrintlnMV(GA3E v) 

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	printf("GA3E_PrintlnMV Demo\nr = ");  GA3E_PrintlnMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void GA3E_PrintMV(GA3E v) 

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	printf("GA3E_PrintMV Demo\nr = ");  GA3E_PrintMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_OverBar(GA3E a) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_OverBar(r);

	printf("GA3E_OverBar Demo\n");
	printf("r =              ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_OverBar(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =       ");  GA3E_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_UnderBar(GA3E a) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_UnderBar(r);

	printf("GA3E_UnderBar Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_UnderBar(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	t = GA3E_OverBar(GA3E_UnderBar(r));
	printf("t = GA3E_OverBar(GA3E_UnderBar(r)) =        ");  GA3E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Reverse(GA3E w) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_Reverse(r);

	printf("GA3E_Reverse Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_Reverse(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	t = GA3E_Reverse(GA3E_Reverse(r));
	printf("t = GA3E_Reverse(GA3E_Reverse(r)) =        ");  GA3E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Parity(GA3E w) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_Parity(r);

	printf("GA3E_Parity Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_Parity(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	t = GA3E_Parity(GA3E_Parity(r));
	printf("t = GA3E_Parity(GA3E_Parity(r)) =        ");  GA3E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Transpose(GA3E w) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_Transpose(r);

	printf("GA3E_Transpose Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_Transpose(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	t = GA3E_Transpose(GA3E_Transpose(r));
	printf("t = GA3E_Transpose(GA3E_Transpose(r)) =        ");  GA3E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Conjugation(GA3E w) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_Conjugation(r);

	printf("GA3E_Conjugation Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_Conjugation(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	t = GA3E_Conjugation(GA3E_Conjugation(r));
	printf("t = GA3E_Conjugation(GA3E_Conjugation(r)) =        ");  GA3E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_CliffordConjugation(GA3E w) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_CliffordConjugation(r);

	printf("GA3E_CliffordConjugation Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_CliffordConjugation(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	t = GA3E_CliffordConjugation(GA3E_CliffordConjugation(r));
	printf("t = GA3E_CliffordConjugation(GA3E_CliffordConjugation(r)) =        ");  GA3E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Dual(GA3E w) ; 

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_Dual(r);

	printf("GA3E_Dual Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_Dual(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	t = GA3E_Dual(GA3E_Dual(r));
	printf("t = GA3E_Dual(GA3E_Dual(r)) =        ");  GA3E_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_DorstDual(GA3E a) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_DorstDual(r);

	printf("GA3E_DorstDual Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_Reverse(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_DorstUnDual(GA3E a) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3E_DorstUnDual(r);

	printf("GA3E_DorstUnDual Demo\n");
	printf("r =               ");  GA3E_PrintMV(r);  printf("\n");
	printf("u = GA3E_DorstUnDual(r) = ");  GA3E_PrintMV(u);  printf("\n");
	det = GA3E_Determinant(u);
	if(det == det_ref) printf("GA3E_Determinant is conserved\n"); else printf("GA3E_Determinant is not conserved\n");
	s = GA3E_Product(r,u);
	printf("s = r*u =        ");  GA3E_PrintMV(s);  printf("\n");

	s = GA3E_DorstDual(r);
	u = GA3E_DorstUnDual(s);

	printf("GA3E_DorstUnDual(GA3E_DorstDual(r)) = "); GA3E_PrintMV(u); printf("\n");

	s = GA3E_DorstUnDual(r);
	u = GA3E_DorstDual(s);

	printf("GA3E_DorstDual(GA3E_DorstUnDual(r)) = "); GA3E_PrintMV(u); printf("\n\n");


//////////////////////////////////////////////////////

// GA3E GA3E_Add(GA3E u,GA3E v)

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	u = GA3E_Add(r, s);

	printf("  r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("  s = ");  GA3E_PrintMV(s);  printf("\n");
	printf("r+s = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Subtract(GA3E u, GA3E v)

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	u = GA3E_Subtract(r, s);

	printf("  r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("  s = ");  GA3E_PrintMV(s);  printf("\n");
	printf("r-s = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int GA3E_Equal(GA3E u, GA3E v) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ; 

	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	printf("t = ");  GA3E_PrintMV(t);  printf("\n");
	if (GA3E_Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (GA3E_Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int GA3E_Not_Equal(GA3E u, GA3E v) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ; 

	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	printf("t = ");  GA3E_PrintMV(t);  printf("\n");
	if (GA3E_Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (GA3E_Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Product(GA3E a, GA3E b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Product(r,s);
	printf("u = r*s = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Divide_By_Constant(GA3E u, double a) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	u = GA3E_Divide_By_Constant(r,2);
	printf("u = r/2 = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Wedge(GA3E a, GA3E b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Wedge(r,s);
	printf("u = GA3E_Wedge(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_AntiWedge(GA3E a, GA3E b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_AntiWedge(r,s);
	printf("u = GA3E_AntiWedge(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Regressive(GA3E a, GA3E b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Regressive(r,s);
	printf("u = GA3E_Regressive(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_RegressiveViaFormula(GA3E a, GA3E b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_RegressiveViaFormula(r,s);
	printf("u = GA3E_RegressiveViaFormula(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_LowerRightViaFormula(GA3E a, GA3E b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_LowerRightViaFormula(r,s);
	printf("u = GA3E_LowerRightViaFormula(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Expander(const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Expander(r,s);
	printf("u = GA3E_Expander(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Conserver(const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Conserver(r,s);
	printf("u = GA3E_Conserver(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Shrinker(const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Shrinker(r,s);
	printf("u = GA3E_Shrinker(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// GA3E GA3E_Symmetric(const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Symmetric(r,s);
	printf("u = GA3E_Symmetric(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_AntiSymmetric(const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_AntiSymmetric(r,s);
	printf("u = GA3E_AntiSymmetric(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Inner(const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_Inner(r,s);
	printf("u = GA3E_Inner(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_LeftContraction (const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_LeftContraction(r,s);
	printf("u = GA3E_LeftContraction(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_RightContraction (const GA3E &a, const GA3E &b) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	printf("s = ");  GA3E_PrintMV(s);  printf("\n");
	u = GA3E_RightContraction(r,s);
	printf("u = GA3E_RightContraction(r, s) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex GA3E_Determinant(GA3E A) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	det = GA3E_Determinant(r);
	printf("u = GA3E_Determinant(r) = %10.3e \n\n",det); 


//////////////////////////////////////////////////////

// GA3E GA3E_Adjugate(GA3E V) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	u = GA3E_Adjugate(r);
	printf("u = GA3E_Adjugate(r) = ");  GA3E_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3E GA3E_Reciprocal(GA3E a) ;

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3E_PrintMV(r);  printf("\n");
	u = GA3E_Reciprocal(r);
	printf("u = GA3E_Reciprocal(r) = ");  GA3E_PrintMV(u);  printf("\n\n");
	s = GA3E_Product(r,u);
	printf("s = r*(1/r) = ");  GA3E_PrintMV(s);  printf("\n\n");

/////////////////////////////////  

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	det_ref = GA3E_Determinant(r);
	int i;

	printf("Verify GA3E_Magic transform preserves determinant\n");
	for(i=0; i<6; i++) {
		s = GA3E_Magic(r,i);
		det = GA3E_Determinant(s);
		if(det != det_ref) printf("GA3E_Magic failure at index %d \n", i);
			else printf("+");
	}
	printf("\n\n");

/////////////////////////////////  

	printf("Verify GA3E_Comp transform preserves determinant\n");
	for(i=0; i<32; i++) {
		s = GA3E_Comp(r,i);
		det = GA3E_Determinant(s);
		if(det != det_ref) printf("GA3E_Comp failure at index %d \n", i);
			else printf("+");
	}
	printf("\n\n");

/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	r = GA3E_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3E_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = GA3E_Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

*/





