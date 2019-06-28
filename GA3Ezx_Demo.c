#include <stdio.h>
#include <stdlib.h>

#include "GA3Ezx_Routines.c"

/////////////////////////////////// Demonstrate GA3Ezx Routines ////////////////////////////////


int main(void)
{
	

	GA3Ezx r,s,t,u;

//	GA3Ezx X, Y, Z;

//	X = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
//	Y = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
//	Z = GA3Ezx_Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
//	t = GA3Ezx_Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

	double det_ref, det;
	det_ref = GA3Ezx_Determinant(r);
	printf("\nGA3Ezx Routines Demo\n\n");

//////////////////////////////////////////////////////

// void GA3Ezx_PrintlnMV(GA3Ezx v) 

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	printf("GA3Ezx_PrintlnMV Demo\nr = ");  GA3Ezx_PrintlnMV(r);  printf("\n");

//////////////////////////////////////////////////////

// void GA3Ezx_PrintMV(GA3Ezx v) 

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	printf("GA3Ezx_PrintMV Demo\nr = ");  GA3Ezx_PrintMV(r);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_OverBar(GA3Ezx a) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_OverBar(r);

	printf("GA3Ezx_OverBar Demo\n");
	printf("r =              ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_OverBar(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =       ");  GA3Ezx_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_UnderBar(GA3Ezx a) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_UnderBar(r);

	printf("GA3Ezx_UnderBar Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_UnderBar(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	t = GA3Ezx_OverBar(GA3Ezx_UnderBar(r));
	printf("t = GA3Ezx_OverBar(GA3Ezx_UnderBar(r)) =        ");  GA3Ezx_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Reverse(GA3Ezx w) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_Reverse(r);

	printf("GA3Ezx_Reverse Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_Reverse(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	t = GA3Ezx_Reverse(GA3Ezx_Reverse(r));
	printf("t = GA3Ezx_Reverse(GA3Ezx_Reverse(r)) =        ");  GA3Ezx_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Parity(GA3Ezx w) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_Parity(r);

	printf("GA3Ezx_Parity Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_Parity(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	t = GA3Ezx_Parity(GA3Ezx_Parity(r));
	printf("t = GA3Ezx_Parity(GA3Ezx_Parity(r)) =        ");  GA3Ezx_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Transpose(GA3Ezx w) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_Transpose(r);

	printf("GA3Ezx_Transpose Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_Transpose(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	t = GA3Ezx_Transpose(GA3Ezx_Transpose(r));
	printf("t = GA3Ezx_Transpose(GA3Ezx_Transpose(r)) =        ");  GA3Ezx_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Conjugation(GA3Ezx w) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_Conjugation(r);

	printf("GA3Ezx_Conjugation Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_Conjugation(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	t = GA3Ezx_Conjugation(GA3Ezx_Conjugation(r));
	printf("t = GA3Ezx_Conjugation(GA3Ezx_Conjugation(r)) =        ");  GA3Ezx_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_CliffordConjugation(GA3Ezx w) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_CliffordConjugation(r);

	printf("GA3Ezx_CliffordConjugation Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_CliffordConjugation(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	t = GA3Ezx_CliffordConjugation(GA3Ezx_CliffordConjugation(r));
	printf("t = GA3Ezx_CliffordConjugation(GA3Ezx_CliffordConjugation(r)) =        ");  GA3Ezx_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Dual(GA3Ezx w) ; 

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_Dual(r);

	printf("GA3Ezx_Dual Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_Dual(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	t = GA3Ezx_Dual(GA3Ezx_Dual(r));
	printf("t = GA3Ezx_Dual(GA3Ezx_Dual(r)) =        ");  GA3Ezx_PrintMV(t);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_DorstDual(GA3Ezx a) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_DorstDual(r);

	printf("GA3Ezx_DorstDual Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_Reverse(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_DorstUnDual(GA3Ezx a) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  

	u = GA3Ezx_DorstUnDual(r);

	printf("GA3Ezx_DorstUnDual Demo\n");
	printf("r =               ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("u = GA3Ezx_DorstUnDual(r) = ");  GA3Ezx_PrintMV(u);  printf("\n");
	det = GA3Ezx_Determinant(u);
	if(det == det_ref) printf("GA3Ezx_Determinant is conserved\n"); else printf("GA3Ezx_Determinant is not conserved\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*u =        ");  GA3Ezx_PrintMV(s);  printf("\n");

	s = GA3Ezx_DorstDual(r);
	u = GA3Ezx_DorstUnDual(s);

	printf("GA3Ezx_DorstUnDual(GA3Ezx_DorstDual(r)) = "); GA3Ezx_PrintMV(u); printf("\n");

	s = GA3Ezx_DorstUnDual(r);
	u = GA3Ezx_DorstDual(s);

	printf("GA3Ezx_DorstDual(GA3Ezx_DorstUnDual(r)) = "); GA3Ezx_PrintMV(u); printf("\n\n");


//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Add(GA3Ezx u,GA3Ezx v)

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	u = GA3Ezx_Add(r, s);

	printf("  r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("  s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	printf("r+s = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Subtract(GA3Ezx u, GA3Ezx v)

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	u = GA3Ezx_Subtract(r, s);

	printf("  r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("  s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	printf("r-s = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// int GA3Ezx_Equal(GA3Ezx u, GA3Ezx v) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ; 

	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	printf("t = ");  GA3Ezx_PrintMV(t);  printf("\n");
	if (GA3Ezx_Equal(r, s)) printf("r == s\n");   else printf("r != s (as expected)\n");
	if (GA3Ezx_Equal(s, t)) printf("s == t (as expected)\n");   else printf("s != t \n");
	printf("\n");

//////////////////////////////////////////////////////

// int GA3Ezx_Not_Equal(GA3Ezx u, GA3Ezx v) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ; 

	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	printf("t = ");  GA3Ezx_PrintMV(t);  printf("\n");
	if (GA3Ezx_Not_Equal(r, s)) printf("r != s (as expected)\n");   else printf("r == s \n");
	if (GA3Ezx_Not_Equal(s, t)) printf("s != t \n");   else printf("s == t (as expected)\n");
	printf("\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Product(GA3Ezx a, GA3Ezx b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Product(r,s);
	printf("u = r*s = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Divide_By_Constant(GA3Ezx u, double a) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	u = GA3Ezx_Divide_By_Constant(r,2);
	printf("u = r/2 = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Wedge(GA3Ezx a, GA3Ezx b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Wedge(r,s);
	printf("u = GA3Ezx_Wedge(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_AntiWedge(GA3Ezx a, GA3Ezx b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_AntiWedge(r,s);
	printf("u = GA3Ezx_AntiWedge(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Regressive(GA3Ezx a, GA3Ezx b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Regressive(r,s);
	printf("u = GA3Ezx_Regressive(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_RegressiveViaFormula(GA3Ezx a, GA3Ezx b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_RegressiveViaFormula(r,s);
	printf("u = GA3Ezx_RegressiveViaFormula(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_LowerRightViaFormula(GA3Ezx a, GA3Ezx b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_LowerRightViaFormula(r,s);
	printf("u = GA3Ezx_LowerRightViaFormula(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Expander(const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Expander(r,s);
	printf("u = GA3Ezx_Expander(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Conserver(const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Conserver(r,s);
	printf("u = GA3Ezx_Conserver(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Shrinker(const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Shrinker(r,s);
	printf("u = GA3Ezx_Shrinker(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

///////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Symmetric(const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Symmetric(r,s);
	printf("u = GA3Ezx_Symmetric(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_AntiSymmetric(const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_AntiSymmetric(r,s);
	printf("u = GA3Ezx_AntiSymmetric(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Inner(const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_Inner(r,s);
	printf("u = GA3Ezx_Inner(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_LeftContraction (const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_LeftContraction(r,s);
	printf("u = GA3Ezx_LeftContraction(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_RightContraction (const GA3Ezx &a, const GA3Ezx &b) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	printf("s = ");  GA3Ezx_PrintMV(s);  printf("\n");
	u = GA3Ezx_RightContraction(r,s);
	printf("u = GA3Ezx_RightContraction(r, s) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// ex GA3Ezx_Determinant(GA3Ezx A) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	det = GA3Ezx_Determinant(r);
	printf("u = GA3Ezx_Determinant(r) = %10.3e \n\n",det); 


//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Adjugate(GA3Ezx V) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	u = GA3Ezx_Adjugate(r);
	printf("u = GA3Ezx_Adjugate(r) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");

//////////////////////////////////////////////////////

// GA3Ezx GA3Ezx_Reciprocal(GA3Ezx a) ;

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	printf("r = ");  GA3Ezx_PrintMV(r);  printf("\n");
	u = GA3Ezx_Reciprocal(r);
	printf("u = GA3Ezx_Reciprocal(r) = ");  GA3Ezx_PrintMV(u);  printf("\n\n");
	s = GA3Ezx_Product(r,u);
	printf("s = r*(1/r) = ");  GA3Ezx_PrintMV(s);  printf("\n\n");

/////////////////////////////////  

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ;  
	det_ref = GA3Ezx_Determinant(r);
	int i;

	printf("Verify GA3Ezx_Magic transform preserves determinant\n");
	for(i=0; i<6; i++) {
		s = GA3Ezx_Magic(r,i);
		det = GA3Ezx_Determinant(s);
		if(det != det_ref) printf("GA3Ezx_Magic failure at index %d \n", i);
			else printf("+");
	}
	printf("\n\n");

/////////////////////////////////  

	printf("Verify GA3Ezx_Comp transform preserves determinant\n");
	for(i=0; i<32; i++) {
		s = GA3Ezx_Comp(r,i);
		det = GA3Ezx_Determinant(s);
		if(det != det_ref) printf("GA3Ezx_Comp failure at index %d \n", i);
			else printf("+");
	}
	printf("\n\n");

/////////////////////////////////  

	return 0;
}

/* copy and paste bin

	r = GA3Ezx_Set(  3,    5,  7, 11,     13, 17, 19,   23) ; 
	s = GA3Ezx_Set( 29,   31, 37, 41,     43, 47, 53,   59) ;
	t = GA3Ezx_Set( 61,   67, 71, 73,     79, 83, 89,   97) ; 

*/





