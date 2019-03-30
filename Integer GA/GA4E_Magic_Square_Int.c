#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "GA4E_Int.c"

///////////////////////////////////

void PrintSimpleBlade(GA4E a) 
{
	int Count = 0;

	if(a.q ==  1) printf(" q     ");	if(a.q == -1) printf("-q     ");	if(a.q ==  0) Count++;

	if(a.x ==  1) printf(" x     ");	if(a.x == -1) printf("-x     ");	if(a.x ==  0) Count++;
	if(a.y ==  1) printf(" y     ");	if(a.y == -1) printf("-y     ");	if(a.y ==  0) Count++;
	if(a.z ==  1) printf(" z     ");	if(a.z == -1) printf("-z     ");	if(a.z ==  0) Count++;
	if(a.t ==  1) printf(" t     ");	if(a.t == -1) printf("-t     ");	if(a.t ==  0) Count++;

	if(a.xy ==  1) printf(" xy    ");	if(a.xy == -1) printf("-xy    ");	if(a.xy ==  0) Count++;
	if(a.xz ==  1) printf(" xz    ");	if(a.xz == -1) printf("-xz    ");	if(a.xz ==  0) Count++;
	if(a.xt ==  1) printf(" xt    ");	if(a.xt == -1) printf("-xt    ");	if(a.xt ==  0) Count++;
	if(a.yz ==  1) printf(" yz    ");	if(a.yz == -1) printf("-yz    ");	if(a.yz ==  0) Count++;
	if(a.yt ==  1) printf(" yt    ");	if(a.yt == -1) printf("-yt    ");	if(a.yt ==  0) Count++;
	if(a.zt ==  1) printf(" zt    ");	if(a.zt == -1) printf("-zt    ");	if(a.zt ==  0) Count++;

	if(a.xyz ==  1) printf(" xyz   ");	if(a.xyz == -1) printf("-xyz   ");	if(a.xyz ==  0) Count++;
	if(a.xyt ==  1) printf(" xyt   ");	if(a.xyt == -1) printf("-xyt   ");	if(a.xyt ==  0) Count++;
	if(a.xzt ==  1) printf(" xzt   ");	if(a.xzt == -1) printf("-xzt   ");	if(a.xzt ==  0) Count++;
	if(a.yzt ==  1) printf(" yzt   ");	if(a.yzt == -1) printf("-yzt   ");	if(a.yzt ==  0) Count++;


	if(a.xyzt ==  1) printf(" xyzt  ");	if(a.xyzt == -1) printf("-xyzt  ");	if(a.xyzt ==  0) Count++;

	if(Count == 16) printf(" 0     ");
	if(Count < 15)  printf("error  ");

}

// ******************************************************

GA4E Multiply_By_Constant(long b, GA4E u)
{
	GA4E a;
	a.q = u.q*b;

	a.x = u.x*b;
	a.y = u.y*b;
	a.z = u.z*b;
	a.t = u.t*b;

	a.xy = u.xy*b;
	a.xz = u.xz*b;
	a.yz = u.yz*b;
	a.xt = u.xt*b;
	a.yt = u.yt*b;
	a.zt = u.zt*b;

	a.xyz = u.xyz*b;
	a.xyt = u.xyt*b;
	a.xzt = u.xzt*b;
	a.yzt = u.yzt*b;

	a.xyzt = u.xyzt*b;

	return a;
}

//////////////////////////////////////////////////////

void printMV(GA4E v) 
{
	printf("(%3ld,  %3ld,%3ld,%3ld,%3ld,  %3ld,%3ld,%3ld, %3ld,%3ld,%3ld,  %3ld,%3ld,%3ld,%3ld, %3ld) ",
		v.q, v.x,v.y,v.z,v.t, v.xy,v.xz,v.yz,v.xt,v.yt,v.zt,  v.xyz,v.xyt,v.xzt,v.yzt,  v.xyzt);
}


///////////////////////////////////

int main(void)
{
	
	int i, ii, jj, kk, mm ;

	time_t start_time, end_time;
	long deltat;

	long a, b,c,d,e,  f,g,h,j,k,l,  m,n,p,r,  s;

	GA4E q, x,y,z,t, xy,xz,xt,yz,yt,zt,  xyz, xyt, xzt, yzt, xyzt;

	GA4E Q, X,Y,Z,T, XY,XZ,XT,YZ,YT,ZT,  XYZ, XYT, XZT, YZT, XYZT;

	time(&start_time);

// initialize the convenience variables.
	
	q = Zero() ; 
	x = Zero() ; 
	y = Zero() ; 
	z = Zero() ; 
	t = Zero() ; 
	xy = Zero() ; 
	xz = Zero() ; 
	xt = Zero() ; 
	yz = Zero() ; 
	yt = Zero() ; 
	zt = Zero() ; 
	xyz = Zero() ; 
	xyt = Zero() ; 
	xzt = Zero() ; 
	yzt = Zero() ; 
	xyzt = Zero() ; 
	
	q.q = 1;	x.x = 1;	y.y = 1;	z.z = 1;	t.t = 1;
	xy.xy = 1;	xz.xz = 1;	xt.xt = 1;	yz.yz = 1;	yt.yt = 1;	zt.zt = 1;
	xyz.xyz = 1;	xyt.xyt = 1;	xzt.xzt = 1;	yzt.yzt = 1;	xyzt.xyzt = 1;

	GA4E Blade[16];

	i = 0;
	Blade[i++] = q;	   Blade[i++] = x;    Blade[i++] = y;    Blade[i++] = z;    Blade[i++] = t;
	Blade[i++] = xy;   Blade[i++] = xz;   Blade[i++] = xt;	 Blade[i++] = yz;   Blade[i++] = yt;   Blade[i++] = zt;
	Blade[i++] = xyz;  Blade[i++] = xyt;  Blade[i++] = xzt;  Blade[i++] = yzt;  Blade[i++] = xyzt;
	
	GA4E Plus[ 6];

	Plus[ 0] =  q      ; 	Plus[ 1] =  x      ; 	Plus[ 2] =  y      ; 	Plus[ 3] =  z      ; 	Plus[ 4] =  t      ; 
	Plus[ 5] =  xyzt   ; 
	
	GA4E MV1,MV2;

	MV1 = Set(  3,    5,  7, 11, 13,   17, 19, 23, 29, 31, 37,    41, 43, 47, 53,  59) ; 

	long det_ref, det;

////////////////////////////////////

//	identify squares of blades

	ii = 0;	jj = 0;
	for(i=0;i<16;i++) {
		MV1 = Blade[i];   
		MV2 = Product(MV1,MV1);
		if (Equal(MV2, q)) {
			printf("Plus[%2d] = ", ii++);
			PrintSimpleBlade(MV1);
			printf(" ; \n");
		} 
	}
	printf("\n");
	ii = 0;	jj = 0;
	for(i=0;i<16;i++) {
		MV1 = Blade[i];   
		MV2 = Product(MV1,MV1);
		if(Not_Equal(MV2, q)) {
			printf("Minus[%2d] = ", jj++);
			PrintSimpleBlade(MV1);
			printf(" ; \n");
		}
	}
	printf("\n");


////////////////////////////////////


	MV1 = Set(  3,    5,  7, 11, 13,   17, 19, 23, 29, 31, 37,    41, 43, 47, 53,  59) ; 
	a = 3; b = 5;  c = 7; d = 11; e = 13;  f = 17;
	g = 19; h = 23; j = 29; k = 31; l = 37;    m = 41; n = 43; p = 47; r =  53; s =  59; 

	det_ref = Determinant(MV1);

// 	build the search loop

	int count = 0;
	for (ii=1;ii< 6;ii++) {
		Q = q;	// always,
		X = Plus[ii];
		for (jj=1;jj< 6;jj++) {
			if(jj != ii) {
				Y = Plus[jj];
				for (kk = 1; kk< 6; kk++) {
					if( (kk != ii) && (kk !=jj) ) {
						Z = Plus[kk];
						for (mm = 1; mm< 6; mm++) {
							if( (mm != ii) && (mm != jj) && (mm != kk) ) {
								T = Plus[mm];
	XY = Product(X,Y) ;	XZ = Product(X,Z) ;	XT = Product(X,T) ;	YZ = Product(Y,Z) ;	YT = Product(Y,T) ;	ZT = Product(Z,T) ;
	XYZ = Product(XY,Z) ;	XYT = Product(XY,T) ;	XZT = Product(XZ,T) ;	YZT = Product(YZ,T) ;	XYZT = Product(XYZ,T) ;

//	MV2 = a*Q + b*X + c*Y + d*Z + e*T
//		+ f*XY  + g*XZ  + h*YZ  + j*XT  + k*YT  + l*ZT
//		+ m*XYZ + n*XYT + p*XZT + r*YZT + s*XYZT;

	MV2 = Multiply_By_Constant(a,Q) ; 

	MV2 = Add( MV2, Multiply_By_Constant(b,X)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(c,Y)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(d,Z)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(e,T)) ; 

	MV2 = Add( MV2, Multiply_By_Constant(f,XY)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(g,XZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(h,YZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(j,XT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(k,YT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(l,ZT)) ; 

	MV2 = Add( MV2, Multiply_By_Constant(m,XYZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(n,XYT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(p,XZT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(r,YZT)) ; 

	MV2 = Add( MV2, Multiply_By_Constant(s,XYZT)) ; 


	det = Determinant(MV2);
	if (det == det_ref)  {
		printf("if(i == %d) MV = ", count++);
		printMV(MV2);
		printf(" ; \n");
		
//		printf("%3d  (",count++);
//		printf("  %10ld  %10ld \n", det, det_ref);
//		PrintSimpleBlade(Q); printf(", ");
//		PrintSimpleBlade(X); printf(", ");
//		PrintSimpleBlade(Y); printf(", ");
//		PrintSimpleBlade(Z); printf(", ");
//		PrintSimpleBlade(T); printf(") ; \n ");
	}

							}
						}
					}
				}
			}
		}
	}


	time(&end_time);
	deltat = difftime(end_time, start_time);
	printf("\nRun time = %ld seconds\n", deltat);


////////////////////////////////////

	return 0;
}

/* copy and paste bin
	a = 3; b = 5;  c = 7; d = 11; e = 13;  f = 17;
	g = 19; h = 23; j = 29; k = 31; l = 37;    m = 41; n = 43; p = 47; r =  53; s =  59; 
	S = 61; R = 67; P = 71; N = 73; M = 79;    L = 83; K = 89; J = 97; H = 101; G = 103; 
	F = 107; E = 109; D = 113; C = 127; B = 131;  A = 137;

	MV1 = GA4E(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA4E(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV1 = GA4E(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	X = GA4E(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 


*/





