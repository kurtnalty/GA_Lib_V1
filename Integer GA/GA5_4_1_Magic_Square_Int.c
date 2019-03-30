#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "GA5_4_1_Int.c"

///////////////////////////////////

void PrintSimpleBlade(GA5_4_1 a) 
{
	int Count = 0;

	if(a.q ==  1) printf(" q     ");	if(a.q == -1) printf("-q     ");	if(a.q ==  0) Count++;

	if(a.w ==  1) printf(" w     ");	if(a.w == -1) printf("-w     ");	if(a.w ==  0) Count++;
	if(a.x ==  1) printf(" x     ");	if(a.x == -1) printf("-x     ");	if(a.x ==  0) Count++;
	if(a.y ==  1) printf(" y     ");	if(a.y == -1) printf("-y     ");	if(a.y ==  0) Count++;
	if(a.z ==  1) printf(" z     ");	if(a.z == -1) printf("-z     ");	if(a.z ==  0) Count++;
	if(a.t ==  1) printf(" t     ");	if(a.t == -1) printf("-t     ");	if(a.t ==  0) Count++;

	if(a.wx ==  1) printf(" wx    ");	if(a.wx == -1) printf("-wx    ");	if(a.wx ==  0) Count++;
	if(a.wy ==  1) printf(" wy    ");	if(a.wy == -1) printf("-wy    ");	if(a.wy ==  0) Count++;
	if(a.wz ==  1) printf(" wz    ");	if(a.wz == -1) printf("-wz    ");	if(a.wz ==  0) Count++;
	if(a.wt ==  1) printf(" wt    ");	if(a.wt == -1) printf("-wt    ");	if(a.wt ==  0) Count++;
	if(a.xy ==  1) printf(" xy    ");	if(a.xy == -1) printf("-xy    ");	if(a.xy ==  0) Count++;
	if(a.xz ==  1) printf(" xz    ");	if(a.xz == -1) printf("-xz    ");	if(a.xz ==  0) Count++;
	if(a.xt ==  1) printf(" xt    ");	if(a.xt == -1) printf("-xt    ");	if(a.xt ==  0) Count++;
	if(a.yz ==  1) printf(" yz    ");	if(a.yz == -1) printf("-yz    ");	if(a.yz ==  0) Count++;
	if(a.yt ==  1) printf(" yt    ");	if(a.yt == -1) printf("-yt    ");	if(a.yt ==  0) Count++;
	if(a.zt ==  1) printf(" zt    ");	if(a.zt == -1) printf("-zt    ");	if(a.zt ==  0) Count++;

	if(a.wxy ==  1) printf(" wxy   ");	if(a.wxy == -1) printf("-wxy   ");	if(a.wxy ==  0) Count++;
	if(a.wxz ==  1) printf(" wxz   ");	if(a.wxz == -1) printf("-wxz   ");	if(a.wxz ==  0) Count++;
	if(a.wxt ==  1) printf(" wxt   ");	if(a.wxt == -1) printf("-wxt   ");	if(a.wxt ==  0) Count++;
	if(a.wyz ==  1) printf(" wyz   ");	if(a.wyz == -1) printf("-wyz   ");	if(a.wyz ==  0) Count++;
	if(a.wyt ==  1) printf(" wyt   ");	if(a.wyt == -1) printf("-wyt   ");	if(a.wyt ==  0) Count++;
	if(a.wzt ==  1) printf(" wzt   ");	if(a.wzt == -1) printf("-wzt   ");	if(a.wzt ==  0) Count++;
	if(a.xyz ==  1) printf(" xyz   ");	if(a.xyz == -1) printf("-xyz   ");	if(a.xyz ==  0) Count++;
	if(a.xyt ==  1) printf(" xyt   ");	if(a.xyt == -1) printf("-xyt   ");	if(a.xyt ==  0) Count++;
	if(a.xzt ==  1) printf(" xzt   ");	if(a.xzt == -1) printf("-xzt   ");	if(a.xzt ==  0) Count++;
	if(a.yzt ==  1) printf(" yzt   ");	if(a.yzt == -1) printf("-yzt   ");	if(a.yzt ==  0) Count++;


	if(a.wxyz ==  1) printf(" wxyz  ");	if(a.wxyz == -1) printf("-xwyz  ");	if(a.wxyz ==  0) Count++;
	if(a.wxyt ==  1) printf(" wxyt  ");	if(a.wxyt == -1) printf("-xwyt  ");	if(a.wxyt ==  0) Count++;
	if(a.wxzt ==  1) printf(" wxzt  ");	if(a.wxzt == -1) printf("-wxzt  ");	if(a.wxzt ==  0) Count++;
	if(a.wyzt ==  1) printf(" wyzt  ");	if(a.wyzt == -1) printf("-wyzt  ");	if(a.wyzt ==  0) Count++;
	if(a.xyzt ==  1) printf(" xyzt  ");	if(a.xyzt == -1) printf("-xyzt  ");	if(a.xyzt ==  0) Count++;

	if(a.wxyzt ==  1) printf(" wxyzt ");	if(a.wxyzt == -1) printf("-wxyzt ");	if(a.wxyzt ==  0) Count++;

	if(Count == 32) printf(" 0     ");
	if(Count < 31)  printf("error  ");

}

// ******************************************************

GA5_4_1 Multiply_By_Constant(long b, GA5_4_1 u)
{
	GA5_4_1 a;
	a.q = u.q*b;

	a.w = u.w*b;
	a.x = u.x*b;
	a.y = u.y*b;
	a.z = u.z*b;
	a.t = u.t*b;

	a.wx = u.wx*b;
	a.wy = u.wy*b;
	a.wz = u.wz*b;
	a.wt = u.wt*b;
	a.xy = u.xy*b;
	a.xz = u.xz*b;
	a.yz = u.yz*b;
	a.xt = u.xt*b;
	a.yt = u.yt*b;
	a.zt = u.zt*b;

	a.wxy = u.wxy*b;
	a.wxz = u.wxz*b;
	a.wxt = u.wxt*b;
	a.wyz = u.wyz*b;
	a.wyt = u.wyt*b;
	a.wzt = u.wzt*b;
	a.xyz = u.xyz*b;
	a.xyt = u.xyt*b;
	a.xzt = u.xzt*b;
	a.yzt = u.yzt*b;

	a.wxyz = u.wxyz*b;
	a.wxyt = u.wxyt*b;
	a.wxzt = u.wxzt*b;
	a.wyzt = u.wyzt*b;
	a.xyzt = u.xyzt*b;

	a.wxyzt = u.wxyzt*b;

	return a;
}

///////////////////////////////////

int main(void)
{
	
	int i, ii, jj, kk, mm, nn ;

	time_t start_time, end_time;
	long deltat;

	long a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A;

	GA5_4_1 q, w,x,y,z,t, wx,wy,wz,wt,xy,xz,xt,yz,yt,zt;
	GA5_4_1 wxy, wxz, wxt, wyz, wyt, wzt, xyz, xyt, xzt, yzt;
	GA5_4_1 wxyz, wxyt, wxzt, wyzt, xyzt, wxyzt;

	GA5_4_1 Q, W,X,Y,Z,T, WX,WY,WZ,WT,XY,XZ,XT,YZ,YT,ZT;
	GA5_4_1 WXY, WXZ, WXT, WYZ, WYT, WZT, XYZ, XYT, XZT, YZT;
	GA5_4_1 WXYZ, WXYT, WXZT, WYZT, XYZT, WXYZT;

	time(&start_time);

// initialize the convenience variables.
	
	q.q = 1;	w.w = 1;	x.x = 1;	y.y = 1;	z.z = 1;	t.t = 1;
	wx.wx = 1;	wy.wy = 1;	wz.wz = 1;	wt.wt = 1;	xy.xy = 1;
	xz.xz = 1;	xt.xt = 1;	yz.yz = 1;	yt.yt = 1;	zt.zt = 1;
	wxy.wxy = 1;	wxz.wxz = 1;	wxt.wxt = 1;	wyz.wyz = 1;	wyt.wyt = 1;
	wzt.wzt = 1;	xyz.xyz = 1;	xyt.xyt = 1;	xzt.xzt = 1;	yzt.yzt = 1;
	wxyz.wxyz = 1;	wxyt.wxyt = 1;	wxzt.wxzt = 1;	wyzt.wyzt = 1;	xyzt.xyzt = 1;	wxyzt.wxyzt = 1;
	GA5_4_1 Blade[32];

	i = 0;
	Blade[i++] = q;	   Blade[i++] = w;    Blade[i++] = x;    Blade[i++] = y;    Blade[i++] = z;   Blade[i++] = t;
	Blade[i++] = wx;   Blade[i++] = wy;   Blade[i++] = wz;   Blade[i++] = wt;   Blade[i++] = xy;
	Blade[i++] = xz;   Blade[i++] = xt;   Blade[i++] = yz;   Blade[i++] = yt;   Blade[i++] = zt;
	Blade[i++] = wxy;  Blade[i++] = wxz;  Blade[i++] = wxt;  Blade[i++] = wyz;  Blade[i++] = wyt;
	Blade[i++] = wzt;  Blade[i++] = xyz;  Blade[i++] = xyt;  Blade[i++] = xzt;  Blade[i++] = yzt;
	Blade[i++] = wxyz; Blade[i++] = wxyt; Blade[i++] = wxzt; Blade[i++] = wyzt; Blade[i++] = xyzt; Blade[i++] = wxyzt;
	
	GA5_4_1 Plus[16], Minus[16];

	Plus[ 0] =  q      ; 	Plus[ 1] =  w      ; 	Plus[ 2] =  x      ; 	Plus[ 3] =  y      ; 	Plus[ 4] =  z      ; 
	Plus[ 5] =  wt     ; 	Plus[ 6] =  xt     ; 	Plus[ 7] =  yt     ; 	Plus[ 8] =  zt     ; 	Plus[ 9] =  wxt    ; 
	Plus[10] =  wyt    ; 	Plus[11] =  wzt    ; 	Plus[12] =  xyt    ; 	Plus[13] =  xzt    ; 	Plus[14] =  yzt    ; 
	Plus[15] =  wxyz   ; 
	
	Minus[ 0] =  t      ; 	Minus[ 1] =  wx     ; 	Minus[ 2] =  wy     ; 	Minus[ 3] =  wz     ; 	Minus[ 4] =  xy     ; 
	Minus[ 5] =  xz     ; 	Minus[ 6] =  yz     ; 	Minus[ 7] =  wxy    ; 	Minus[ 8] =  wxz    ; 	Minus[ 9] =  wyz    ; 
	Minus[10] =  xyz    ; 	Minus[11] =  wxyt   ; 	Minus[12] =  wxzt   ; 	Minus[13] =  wyzt   ; 	Minus[14] =  xyzt   ; 
	Minus[15] =  wxyzt  ; 
	

	GA5_4_1 MV1,MV2;


	MV1 = Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	a = 3; b = 5;  c = 7; d = 11; e = 13;  f = 17;
	g = 19; h = 23; j = 29; k = 31; l = 37;    m = 41; n = 43; p = 47; r =  53; s =  59; 
	S = 61; R = 67; P = 71; N = 73; M = 79;    L = 83; K = 89; J = 97; H = 101; G = 103; 
	F = 107; E = 109; D = 113; C = 127; B = 131;  A = 137;
	complex det_ref, det;

////////////////////////////////////

//	identify squares of blades

	ii = 0;	jj = 0;
	for(i=0;i<32;i++) {
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
	for(i=0;i<32;i++) {
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

	MV1 = Set(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 

	det_ref = Determinant(MV1);

// 	build the search loop

	int count = 0;
	for (ii=1;ii<16;ii++) {
		Q = q;	// always,
		W = Plus[ii];
		for (jj=1;jj<16;jj++) {
			if(jj != ii) {
				X = Plus[jj];
				for (kk = 1; kk<16; kk++) {
					if( (kk != ii) && (kk !=jj) ) {
						Y = Plus[kk];
						for (mm = 1; mm<16; mm++) {
							if( (mm != ii) && (mm != jj) && (mm != kk) ) {
								Z = Plus[mm];
								for (nn = 0; nn < 16; nn++) {
									T = Minus[nn];
	WX = Product(W,X) ;	WY = Product(W,Y) ;	WZ = Product(W,Z) ;	WT = Product(W,T) ;	XY = Product(X,Y) ;
	XZ = Product(X,Z) ;	XT = Product(X,T) ;	YZ = Product(Y,Z) ;	YT = Product(Y,T) ;	ZT = Product(Z,T) ;
	WXY = Product(WX,Y) ;	WXZ = Product(WX,Z) ;	WXT = Product(WX,T) ;	WYZ = Product(WY,Z) ;	WYT = Product(WY,T) ;
	WZT = Product(WZ,T) ;	XYZ = Product(XY,Z) ;	XYT = Product(XY,T) ;	XZT = Product(XZ,T) ;	YZT = Product(YZ,T) ;
	WXYZ = Product(WXY,Z) ;	WXYT = Product(WXY,T) ;	WXZT = Product(WXZ,T) ;	WYZT = Product(WYZ,T) ;	XYZT = Product(XYZ,T) ;
	WXYZT = Product(WXYZ,T) ;

//	MV2 = a*Q + b*W + c*X + d*Y + e*Z + f*T
//		+ g*WX  + h*WY  + j*WZ  + k*WT  + l*XY  + m*XZ  + n*XT  + p*YZ  + r*YT  + s*ZT
//		+ S*WXY + R*WXZ + P*WXT + N*WYZ + M*WYT + L*WZT + K*XYZ + J*XYT + H*XZT + G*YZT
//		+ F*WXYZ + E*WXYT + D*WXZT + C*WYZT + B*XYZT + A*WXYZT;

	MV2 = Multiply_By_Constant(a,Q) ; 

	MV2 = Add( MV2, Multiply_By_Constant(b,W)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(c,X)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(d,Y)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(e,Z)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(f,T)) ; 

	MV2 = Add( MV2, Multiply_By_Constant(g,WX)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(h,WY)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(j,WZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(k,WT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(l,XY)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(m,XZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(n,XT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(p,YZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(r,YT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(s,ZT)) ; 

	MV2 = Add( MV2, Multiply_By_Constant(S,WXY)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(R,WXZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(P,WXT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(N,WYZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(M,WYT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(L,WZT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(K,XYZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(J,XYT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(H,XZT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(G,YZT)) ; 

	MV2 = Add( MV2, Multiply_By_Constant(F,WXYZ)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(E,WXYT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(D,WXZT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(C,WYZT)) ; 
	MV2 = Add( MV2, Multiply_By_Constant(B,XYZT)) ; 

	MV2 = Add( MV2, Multiply_By_Constant(A,WXYZT)) ; 


	det = Determinant(MV2);
	if((det.r == det_ref.r) && (det.i == det_ref.i)) {
//		printf("if(i == %d) MV = ", count++);
//		PrintMV(MV2);
//		printf(" ; \n");
		
		printf("%3d  (",count++);
		PrintSimpleBlade(Q); printf(", ");
		PrintSimpleBlade(W); printf(", ");
		PrintSimpleBlade(X); printf(", ");
		PrintSimpleBlade(Y); printf(", ");
		PrintSimpleBlade(Z); printf(", ");
		PrintSimpleBlade(T); printf(") ; \n ");


	}
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

	MV1 = GA5_4_1(a_q, a_w,a_x,a_y,a_z,a_t, 
		a_wx, a_wy, a_wz, a_wt, a_xy, a_xz, a_xt, a_yz, a_yt, a_zt, 
		a_wxy, a_wxz, a_wxt, a_wyz, a_wyt, a_wzt, a_xyz, a_xyt, a_xzt, a_yzt,
		a_wxyz, a_wxyt, a_wxzt, a_wyzt, a_xyzt, a_wxyzt);

	MV2 = GA5_4_1(b_q, b_w,b_x,b_y,b_z,b_t, 
		b_wx, b_wy, b_wz, b_wt, b_xy, b_xz, b_xt, b_yz, b_yt, b_zt, 
		b_wxy, b_wxz, b_wxt, b_wyz, b_wyt, b_wzt, b_xyz, b_xyt, b_xzt, b_yzt,
		b_wxyz, b_wxyt, b_wxzt, b_wyzt, b_xyzt, b_wxyzt);

	MV1 = GA5_4_1(a, b,c,d,e,f, g,h,j,k,l, m,n,p,r,s, S,R,P,N,M, L,K,J,H,G, F,E,D,C,B, A);

	X = GA5_4_1(  3,    5,  7, 11, 13, 17, 
		 19, 23, 29, 31, 37,    41, 43, 47, 53, 59,
		 61, 67, 71, 73, 79,    83, 89, 97,101,103,
		107,109,113,127,131,   137) ; 



*/





