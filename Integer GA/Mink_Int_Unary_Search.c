#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Mink_Int.c"

/////////////////////////////////// 

Mink Mink_Unary(Mink a, int i)
{
	Mink b;
	int SignArray[16];
	int j,k;

	k = 1;
	for (j=0;j<16;j++) {
		if ( (k & i) != 0) SignArray[j] = -1; else SignArray[j] = +1;
		k <<= 1;
	}

	b.q = SignArray[ 0]*a.q;

	b.x = SignArray[ 1]*a.x;
	b.y = SignArray[ 2]*a.y;
	b.z = SignArray[ 3]*a.z;
	b.t = SignArray[ 4]*a.t;

	b.xy = SignArray[ 5]*a.xy;
	b.xz = SignArray[ 6]*a.xz;
	b.yz = SignArray[ 7]*a.yz;
	b.xt = SignArray[ 8]*a.xt;
	b.yt = SignArray[ 9]*a.yt;
	b.zt = SignArray[10]*a.zt;

	b.xyz = SignArray[11]*a.xyz;
	b.xyt = SignArray[12]*a.xyt;
	b.xzt = SignArray[13]*a.xzt;
	b.yzt = SignArray[14]*a.yzt;

	b.xyzt = SignArray[15]*a.xyzt;


	return b;

}

/////////////////////////////////// 

	int ZCount(Mink a)
{
	int Count = 0;

	if(a.q    == 0) Count++;

	if(a.x    == 0) Count++;
	if(a.y    == 0) Count++;
	if(a.z    == 0) Count++;
	if(a.t    == 0) Count++;

	if(a.xy   == 0) Count++;
	if(a.xz   == 0) Count++;
	if(a.yz   == 0) Count++;
	if(a.xt   == 0) Count++;
	if(a.yt   == 0) Count++;
	if(a.zt   == 0) Count++;

	if(a.xyz  == 0) Count++;
	if(a.xyt  == 0) Count++;
	if(a.xzt  == 0) Count++;
	if(a.yzt  == 0) Count++;

	if(a.xyzt == 0) Count++;

	return Count;

}

/////////////////////////////////// 




int main(void)
{
	int i,Count;
	long det_ref, det_neg, det;
	Mink r,s,t,u;

	r = Set(  3,    5,  7, 11, 13,     17, 19, 23, 29, 31, 37,    41, 43, 47, 53,    59) ; 

	det_ref = Determinant(r);
	det_neg = -det_ref;

	Count = 0;
	for (i=0;i<0x10000;i++) {
		s = Mink_Unary(r,i);
		det = Determinant(s);
		if (det == det_ref) {
			Count++;
			t = Product(r,s);
			u = Product(s,r);
			printf("dets match at i = %d, Score = %d ",i,ZCount(t));
//			printf("  s = ");  PrintMV(s);  printf("\n");
//			printf("r*s = ");  PrintMV(t);  printf("\n");
//			printf("s*r = ");  PrintMV(u);  printf("\n");
			if(Equal(t,u)) printf("Commutes");
			printf("\n");
		}
	}

	printf("\nCount = %d\n\n",Count);


	Count = 0;
	for (i=0;i<0x10000;i++) {
		s = Mink_Unary(r,i);
		det = Determinant(s);
		if (det == det_neg) {
			Count++;
			t = Product(r,s);
			u = Product(s,r);
			printf("dets match at i = %d, Score = %d ",i,ZCount(t));
//			printf("  s = ");  PrintMV(s);  printf("\n");
//			printf("r*s = ");  PrintMV(t);  printf("\n");
//			printf("s*r = ");  PrintMV(u);  printf("\n");
			if(Equal(t,u)) printf("Commutes");
			printf("\n");
		}
	}

	printf("\nCount = %d\n\n",Count);

//////////////////////////////////////////////////////

	return 0;
}






