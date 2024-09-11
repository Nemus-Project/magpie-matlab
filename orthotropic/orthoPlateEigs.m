function [Q,Om,Nx,Ny,biHarm] = orthoPlateEigs(rho,Ex,Ey,Gxy,nux,Lx,Ly,Lz,hx,hy,K0y,R0y,Kx0,Rx0,KLy,RLy,KxL,RxL,Nmodes,printfreqs)


%--- derived parameters (don't change here)
nuy     = Ey/Ex*nux ;
D1      = Ex * Lz^3 / 12 / (1-nux*nuy) ;
D3      = Ey * Lz^3 / 12 / (1-nux*nuy) ;
D4      = Gxy * Lz^3 / 3 ;
D2      = nuy*Ex*Lz^3/ 6 / (1-nux*nuy) ;
Nx      = floor(Lx/hx) ;
Ny      = floor(Ly/hy) ;
%----------------------------


% MATRIX BUILDER

%--- build matrix in blocks

a0 = ones(Ny+1,1) ;
a1 = ones(Ny,1) ;
a2 = ones(Ny-1,1) ;



[D00u00,D00u10,D00u20,D00u01,D00u02,D00u11] = D00_coeffs(K0y,R0y,Kx0,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); % D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[D01u01,D01u11,D01u21,D01u00,D01u02,D01u03,D01u12,D01u10] = D01_coeffs(K0y,R0y,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); % D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[D02u02,D02u12,D02u22,D02u01,D02u03,D02u04,D02u00,D02u13,D02u11] = D02_coeffs(K0y,R0y,hx,hy,nuy,D1,D2,D3,D4); % D02_coeffs(K0y,R0y,h,D,nu) ;
[D0Nu0N,D0Nu1N,D0Nu2N,D0Nu0Nm1,D0Nu0Nm2,D0Nu1Nm1] = D00_coeffs(K0y,R0y,KxL,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %v D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,D0Nm1u1Nm1,D0Nm1u2Nm1,D0Nm1u0N,D0Nm1u0Nm2,D0Nm1u0Nm3,D0Nm1u1Nm2,D0Nm1u1N] = D01_coeffs(K0y,R0y,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

%-- Blk11
W0 = D02u02*a0 ; W1 = D02u03*a1; W2 = D02u04*a2 ; Wm1 = D02u01*a1; Wm2 = D02u00*a2 ;

Blk11               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1) + diag(W2,2) + diag(Wm2,-2)) ;
Blk11(1,1)          = D00u00 ;
Blk11(1,2)          = D00u01 ;
Blk11(1,3)          = D00u02 ;
Blk11(2,1)          = D01u00 ;
Blk11(2,2)          = D01u01 ;
Blk11(2,3)          = D01u02 ;
Blk11(2,4)          = D01u03 ;
Blk11(end,end)      = D0Nu0N ;
Blk11(end,end-1)    = D0Nu0Nm1 ;
Blk11(end,end-2)    = D0Nu0Nm2 ;
Blk11(end-1,end)    = D0Nm1u0N ;
Blk11(end-1,end-1)  = D0Nm1u0Nm1 ;
Blk11(end-1,end-2)  = D0Nm1u0Nm2 ;
Blk11(end-1,end-3)  = D0Nm1u0Nm3 ;

%-- Blk12
W0 = D02u12*a0 ; W1 = D02u13*a1 ; Wm1 = D02u11*a1 ;

Blk12               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1))  ;
Blk12(1,1)          = D00u10 ;
Blk12(1,2)          = D00u11 ;
Blk12(2,1)          = D01u10 ;
Blk12(2,2)          = D01u11 ;
Blk12(2,3)          = D01u12 ;
Blk12(end,end)      = D0Nu1N ;
Blk12(end,end-1)    = D0Nu1Nm1 ;
Blk12(end-1,end)    = D0Nm1u1N ;
Blk12(end-1,end-1)  = D0Nm1u1Nm1 ;
Blk12(end-1,end-2)  = D0Nm1u1Nm2 ;


%-- Blk13
W0 = D02u22*a0 ;

Blk13               = sparse(diag(W0))   ;
Blk13(1,1)          = D00u20 ;
Blk13(2,2)          = D01u21 ;
Blk13(end,end)      = D0Nu2N ;
Blk13(end-1,end-1)  = D0Nm1u2Nm1 ;




[D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01] = D10_coeffs(R0y,Kx0,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); % D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[D11u11,D11u12,D11u13,D11u10,D11u01,D11u21,D11u31,D11u22,D11u20,D11u00,D11u02]  = D11_coeffs(R0y,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); %D11_coeffs(R0y,Rx0,h,D,nu) ;
[D12u12,D12u13,D12u14,D12u11,D12u10,D12u02,D12u22,D12u32,D12u23,D12u21,D12u01,D12u03]  = D12_coeffs(R0y,hx,hy,nuy,D1,D2,D3,D4); %D12_coeffs(R0y,h,D,nu) ;
[D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(R0y,KxL,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[D1Nm1u1Nm1,D1Nm1u1Nm2,D1Nm1u1Nm3,D1Nm1u1N,D1Nm1u0Nm1,D1Nm1u2Nm1,D1Nm1u3Nm1,D1Nm1u2Nm2,D1Nm1u2N,D1Nm1u0N,D1Nm1u0Nm2]  = D11_coeffs(R0y,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %D11_coeffs(R0y,RxL,h,D,nu) ;


%-- Blk21
W0 = D12u02*a0 ; W1 = D12u03*a1 ; Wm1 = D12u01*a1 ;

Blk21               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1))  ;
Blk21(1,1)          = D10u00 ;
Blk21(1,2)          = D10u01 ;
Blk21(2,1)          = D11u00 ;
Blk21(2,2)          = D11u01 ;
Blk21(2,3)          = D11u02 ;
Blk21(end,end)      = D1Nu0N ;
Blk21(end,end-1)    = D1Nu0Nm1 ;
Blk21(end-1,end)    = D1Nm1u0N ;
Blk21(end-1,end-1)  = D1Nm1u0Nm1 ;
Blk21(end-1,end-2)  = D1Nm1u0Nm2 ;


%-- Blk22
W0 = D12u12*a0 ; W1 = D12u13*a1 ; W2 = D12u14*a2 ; Wm1 = D12u11*a1 ; Wm2 = D12u10*a2 ;

Blk22               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1) + diag(W2,2) + diag(Wm2,-2)) ;
Blk22(1,1)          = D10u10 ;
Blk22(1,2)          = D10u11 ;
Blk22(1,3)          = D10u12 ;
Blk22(2,1)          = D11u10 ;
Blk22(2,2)          = D11u11 ;
Blk22(2,3)          = D11u12 ;
Blk22(2,4)          = D11u13 ;

Blk22(end,end)      = D1Nu1N ;
Blk22(end,end-1)    = D1Nu1Nm1 ;
Blk22(end,end-2)    = D1Nu1Nm2 ;
Blk22(end-1,end)    = D1Nm1u1N ;
Blk22(end-1,end-1)  = D1Nm1u1Nm1 ;
Blk22(end-1,end-2)  = D1Nm1u1Nm2 ;
Blk22(end-1,end-3)  = D1Nm1u1Nm3 ;


%-- Blk23
W0 = D12u22*a0 ; W1 = D12u23*a1 ; Wm1 = D12u21*a1 ;

Blk23               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1)) ;
Blk23(1,1)          = D10u20 ;
Blk23(1,2)          = D10u21 ;
Blk23(2,1)          = D11u20 ;
Blk23(2,2)          = D11u21 ;
Blk23(2,3)          = D11u22 ;

Blk23(end,end)      = D1Nu2N ;
Blk23(end,end-1)    = D1Nu2Nm1 ;
Blk23(end-1,end)    = D1Nm1u2N ;
Blk23(end-1,end-1)  = D1Nm1u2Nm1 ;
Blk23(end-1,end-2)  = D1Nm1u2Nm2 ;


%-- Blk24
W0 = D12u32*a0 ;

Blk24               = sparse(diag(W0))   ;
Blk24(1,1)          = D10u30 ;
Blk24(2,2)          = D11u31 ;
Blk24(end,end)      = D1Nu3N ;
Blk24(end-1,end-1)  = D1Nm1u3Nm1 ;

[D20u20,D20u21,D20u22,D20u10,D20u30,D20u40,D20u00,D20u31,D20u11] = D20_coeffs(Kx0,Rx0,hx,hy,nux,D1,D2,D3,D4); % D20_coeffs(Kx0,Rx0,h,D,nu) ;
[D21u21,D21u22,D21u23,D21u20,D21u11,D21u31,D21u41,D21u01,D21u32,D21u30,D21u10,D21u12] = D21_coeffs(Rx0,hx,hy,nux,D1,D2,D3,D4); % D21_coeffs(Rx0,h,D,nu) ;
[D22u20,D22u11,D22u21,D22u31,D22u02,D22u12,D22u22,D22u32,D22u42,D22u13,D22u23,D22u33,D22u24] = D22_coeffs(hx,hy,D1,D2,D3,D4) ;
[D2Nu2N,D2Nu2Nm1,D2Nu2Nm2,D2Nu1N,D2Nu3N,D2Nu4N,D2Nu0N,D2Nu3Nm1,D2Nu1Nm1] = D20_coeffs(KxL,RxL,hx,hy,nux,D1,D2,D3,D4); %D20_coeffs(KxL,RxL,h,D,nu) ;
[D2Nm1u2Nm1,D2Nm1u2Nm2,D2Nm1u2Nm3,D2Nm1u2N,D2Nm1u1Nm1,D2Nm1u3Nm1,D2Nm1u4Nm1,D2Nm1u0Nm1,D2Nm1u3Nm2,D2Nm1u3N,D2Nm1u1N,D2Nm1u1Nm2] = D21_coeffs(RxL,hx,hy,nux,D1,D2,D3,D4); %D21_coeffs(RxL,h,D,nu) ;

%-- Blk31
W0 = D22u02*a0 ;

Blk31               = sparse(diag(W0))   ;
Blk31(1,1)          = D20u00 ;
Blk31(2,2)          = D21u01 ;
Blk31(end,end)      = D2Nu0N ;
Blk31(end-1,end-1)  = D2Nm1u0Nm1 ;

%-- Blk32
W0 = D22u12*a0 ; W1 = D22u13*a1 ; Wm1 = D22u11*a1 ;

Blk32               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1)) ;
Blk32(1,1)          = D20u10 ;
Blk32(1,2)          = D20u11 ;
Blk32(2,1)          = D21u10 ;
Blk32(2,2)          = D21u11 ;
Blk32(2,3)          = D21u12 ;

Blk32(end,end)      = D2Nu1N  ;
Blk32(end,end-1)    = D2Nu1Nm1 ;
Blk32(end-1,end)    = D2Nm1u1N ;
Blk32(end-1,end-1)  = D2Nm1u1Nm1 ;
Blk32(end-1,end-2)  = D2Nm1u1Nm2 ;


%-- Blk33
W0 = D22u22*a0 ; W1 = D22u23*a1; W2 = D22u24*a2 ; Wm1 = D22u21*a1; Wm2 = D22u20*a2 ;


Blk33               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1) + diag(W2,2) + diag(Wm2,-2)) ;
Blk33(1,1)          = D20u20 ;
Blk33(1,2)          = D20u21 ;
Blk33(1,3)          = D20u22 ;
Blk33(2,1)          = D21u20 ;
Blk33(2,2)          = D21u21 ;
Blk33(2,3)          = D21u22 ;
Blk33(2,4)          = D21u23 ;

Blk33(end,end)      = D2Nu2N ;
Blk33(end,end-1)    = D2Nu2Nm1 ;
Blk33(end,end-2)    = D2Nu2Nm2 ;
Blk33(end-1,end)    = D2Nm1u2N ;
Blk33(end-1,end-1)  = D2Nm1u2Nm1 ;
Blk33(end-1,end-2)  = D2Nm1u2Nm2 ;
Blk33(end-1,end-3)  = D2Nm1u2Nm3 ;

%-- Blk34
W0 = D22u32*a0 ; W1 = D22u33*a1 ; Wm1 = D22u31*a1 ;

Blk34               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1)) ;
Blk34(1,1)          =  D20u30 ;
Blk34(1,2)          =  D20u31 ;
Blk34(2,1)          =  D21u30 ;
Blk34(2,2)          =  D21u31 ;
Blk34(2,3)          =  D21u32 ;

Blk34(end,end)      =  D2Nu3N;
Blk34(end,end-1)    =  D2Nu3Nm1;
Blk34(end-1,end)    =  D2Nm1u3N;
Blk34(end-1,end-1)  =  D2Nm1u3Nm1;
Blk34(end-1,end-2)  =  D2Nm1u3Nm2;

%-- Blk35
W0 = D22u42*a0 ;

Blk35               = sparse(diag(W0))   ;
Blk35(1,1)          = D20u40 ;
Blk35(2,2)          = D21u41 ;
Blk35(end,end)      = D2Nu4N ;
Blk35(end-1,end-1)  = D2Nm1u4Nm1 ;




[D00u00,D00u10,D00u20,D00u01,D00u02,D00u11] = D00_coeffs(KLy,RLy,Kx0,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); %D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[D01u01,D01u11,D01u21,D01u00,D01u02,D01u03,D01u12,D01u10] = D01_coeffs(KLy,RLy,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); %D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[D02u02,D02u12,D02u22,D02u01,D02u03,D02u04,D02u00,D02u13,D02u11] = D02_coeffs(KLy,RLy,hx,hy,nuy,D1,D2,D3,D4); %D02_coeffs(KLy,RLy,h,D,nu) ;
[D0Nu0N,D0Nu1N,D0Nu2N,D0Nu0Nm1,D0Nu0Nm2,D0Nu1Nm1] = D00_coeffs(KLy,RLy,KxL,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,D0Nm1u1Nm1,D0Nm1u2Nm1,D0Nm1u0N,D0Nm1u0Nm2,D0Nm1u0Nm3,D0Nm1u1Nm2,D0Nm1u1N] = D01_coeffs(KLy,RLy,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

%-- BlkMM
W0 = D02u02*a0 ; W1 = D02u03*a1; W2 = D02u04*a2 ; Wm1 = D02u01*a1; Wm2 = D02u00*a2 ;

BlkMM               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1) + diag(W2,2) + diag(Wm2,-2)) ;
BlkMM(1,1)          = D00u00 ;
BlkMM(1,2)          = D00u01 ;
BlkMM(1,3)          = D00u02 ;
BlkMM(2,1)          = D01u00 ;
BlkMM(2,2)          = D01u01 ;
BlkMM(2,3)          = D01u02 ;
BlkMM(2,4)          = D01u03 ;
BlkMM(end,end)      = D0Nu0N ;
BlkMM(end,end-1)    = D0Nu0Nm1 ;
BlkMM(end,end-2)    = D0Nu0Nm2 ;
BlkMM(end-1,end)    = D0Nm1u0N ;
BlkMM(end-1,end-1)  = D0Nm1u0Nm1 ;
BlkMM(end-1,end-2)  = D0Nm1u0Nm2 ;
BlkMM(end-1,end-3)  = D0Nm1u0Nm3 ;

%-- BlkMMm1
W0 = D02u12*a0 ; W1 = D02u13*a1 ; Wm1 = D02u11*a1 ;

BlkMMm1               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1))  ;
BlkMMm1(1,1)          = D00u10 ;
BlkMMm1(1,2)          = D00u11 ;
BlkMMm1(2,1)          = D01u10 ;
BlkMMm1(2,2)          = D01u11 ;
BlkMMm1(2,3)          = D01u12 ;
BlkMMm1(end,end)      = D0Nu1N ;
BlkMMm1(end,end-1)    = D0Nu1Nm1 ;
BlkMMm1(end-1,end)    = D0Nm1u1N ;
BlkMMm1(end-1,end-1)  = D0Nm1u1Nm1 ;
BlkMMm1(end-1,end-2)  = D0Nm1u1Nm2 ;


%-- BlkMMm2
W0 = D02u22*a0 ;

BlkMMm2               = sparse(diag(W0))   ;
BlkMMm2(1,1)          = D00u20 ;
BlkMMm2(2,2)          = D01u21 ;
BlkMMm2(end,end)      = D0Nu2N ;
BlkMMm2(end-1,end-1)  = D0Nm1u2Nm1 ;



[D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01] = D10_coeffs(RLy,Kx0,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); %D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[D11u11,D11u12,D11u13,D11u10,D11u01,D11u21,D11u31,D11u22,D11u20,D11u00,D11u02] = D11_coeffs(RLy,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4); %D11_coeffs(RLy,Rx0,h,D,nu) ;
[D12u12,D12u13,D12u14,D12u11,D12u10,D12u02,D12u22,D12u32,D12u23,D12u21,D12u01,D12u03] = D12_coeffs(RLy,hx,hy,nuy,D1,D2,D3,D4); %D12_coeffs(RLy,h,D,nu) ;
[D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(RLy,KxL,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[D1Nm1u1Nm1,D1Nm1u1Nm2,D1Nm1u1Nm3,D1Nm1u1N,D1Nm1u0Nm1,D1Nm1u2Nm1,D1Nm1u3Nm1,D1Nm1u2Nm2,D1Nm1u2N,D1Nm1u0N,D1Nm1u0Nm2] = D11_coeffs(RLy,RxL,hx,hy,nux,nuy,D1,D2,D3,D4); %D11_coeffs(RLy,RxL,h,D,nu) ;

%-- BlkMm1M
W0 = D12u02*a0 ; W1 = D12u03*a1 ; Wm1 = D12u01*a1 ;

BlkMm1M               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1))  ;
BlkMm1M(1,1)          = D10u00 ;
BlkMm1M(1,2)          = D10u01 ;
BlkMm1M(2,1)          = D11u00 ;
BlkMm1M(2,2)          = D11u01 ;
BlkMm1M(2,3)          = D11u02 ;
BlkMm1M(end,end)      = D1Nu0N ;
BlkMm1M(end,end-1)    = D1Nu0Nm1 ;
BlkMm1M(end-1,end)    = D1Nm1u0N ;
BlkMm1M(end-1,end-1)  = D1Nm1u0Nm1 ;
BlkMm1M(end-1,end-2)  = D1Nm1u0Nm2 ;


%-- BlkMm1Mm1
W0 = D12u12*a0 ; W1 = D12u13*a1 ; W2 = D12u14*a2 ; Wm1 = D12u11*a1 ; Wm2 = D12u10*a2 ;


BlkMm1Mm1               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1) + diag(W2,2) + diag(Wm2,-2)) ;
BlkMm1Mm1(1,1)          = D10u10 ;
BlkMm1Mm1(1,2)          = D10u11 ;
BlkMm1Mm1(1,3)          = D10u12 ;
BlkMm1Mm1(2,1)          = D11u10 ;
BlkMm1Mm1(2,2)          = D11u11 ;
BlkMm1Mm1(2,3)          = D11u12 ;
BlkMm1Mm1(2,4)          = D11u13 ;

BlkMm1Mm1(end,end)      = D1Nu1N ;
BlkMm1Mm1(end,end-1)    = D1Nu1Nm1 ;
BlkMm1Mm1(end,end-2)    = D1Nu1Nm2 ;
BlkMm1Mm1(end-1,end)    = D1Nm1u1N ;
BlkMm1Mm1(end-1,end-1)  = D1Nm1u1Nm1 ;
BlkMm1Mm1(end-1,end-2)  = D1Nm1u1Nm2 ;
BlkMm1Mm1(end-1,end-3)  = D1Nm1u1Nm3 ;


%-- BlkMm1Mm2
W0 = D12u22*a0 ; W1 = D12u23*a1 ; Wm1 = D12u21*a1 ;

BlkMm1Mm2               = sparse(diag(W0) + diag(W1,1) + diag(Wm1,-1)) ;
BlkMm1Mm2(1,1)          = D10u20 ;
BlkMm1Mm2(1,2)          = D10u21 ;
BlkMm1Mm2(2,1)          = D11u20 ;
BlkMm1Mm2(2,2)          = D11u21 ;
BlkMm1Mm2(2,3)          = D11u22 ;

BlkMm1Mm2(end,end)      = D1Nu2N ;
BlkMm1Mm2(end,end-1)    = D1Nu2Nm1 ;
BlkMm1Mm2(end-1,end)    = D1Nm1u2N ;
BlkMm1Mm2(end-1,end-1)  = D1Nm1u2Nm1 ;
BlkMm1Mm2(end-1,end-2)  = D1Nm1u2Nm2 ;


%-- BlkMm1Mm3
W0 = D12u32*a0 ;

BlkMm1Mm3               = sparse(diag(W0))   ;
BlkMm1Mm3(1,1)          = D10u30 ;
BlkMm1Mm3(2,2)          = D11u31 ;
BlkMm1Mm3(end,end)      = D1Nu3N ;
BlkMm1Mm3(end-1,end-1)  = D1Nm1u3Nm1 ;


biHarm = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;


for m = 3 : Nx - 1
    biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-1)+1 : (Ny+1)*m)        = Blk33 ;
    biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-2)+1 : (Ny+1)*(m-1))    = Blk32 ;
    biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m-3)+1 : (Ny+1)*(m-2))    = Blk31 ;
    biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*m+1 : (Ny+1)*(m+1))        = Blk34 ;
    biHarm((Ny+1)*(m-1)+1 : (Ny+1)*m, (Ny+1)*(m+1)+1 : (Ny+1)*(m+2))    = Blk35 ;
end

biHarm(1:Ny+1,1:Ny+1)           = Blk11 ;
biHarm(1:Ny+1,Ny+2:2*Ny+2)      = Blk12 ;
biHarm(1:Ny+1,2*Ny+3:3*Ny+3)    = Blk13 ;


biHarm(Nx*(Ny+1)+1:(Ny+1)*(Nx+1),Nx*(Ny+1)+1:(Ny+1)*(Nx+1))         = BlkMM ;
biHarm(Nx*(Ny+1)+1:(Ny+1)*(Nx+1),(Nx-1)*(Ny+1)+1:(Ny+1)*Nx)         = BlkMMm1 ;
biHarm(Nx*(Ny+1)+1:(Ny+1)*(Nx+1),(Nx-2)*(Ny+1)+1:(Ny+1)*(Nx-1))     = BlkMMm2 ;


biHarm(Ny+2:2*Ny+2,1:Ny+1)          = Blk21 ;
biHarm(Ny+2:2*Ny+2,Ny+2:2*Ny+2)     = Blk22 ;
biHarm(Ny+2:2*Ny+2,2*Ny+3:3*Ny+3)   = Blk23 ;
biHarm(Ny+2:2*Ny+2,3*Ny+4:4*Ny+4)   = Blk24 ;

biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),Nx*(Ny+1)+1:(Ny+1)*(Nx+1))         = BlkMm1M ;
biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-1)*(Ny+1)+1:(Ny+1)*(Nx))       = BlkMm1Mm1 ;
biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-2)*(Ny+1)+1:(Ny+1)*(Nx-1))     = BlkMm1Mm2 ;
biHarm((Ny+1)*(Nx-1)+1:Nx*(Ny+1),(Nx-3)*(Ny+1)+1:(Ny+1)*(Nx-2))     = BlkMm1Mm3 ;


% EIGENVALUES

[Q,Dm] = eigs(biHarm,Nmodes,'smallestabs') ;
[~,indSort] = sort(diag((Dm))) ;
Q = Q(:,indSort) ;


freqs = sqrt(abs(diag(Dm)))*sqrt(1/rho/Lz)/2/pi ;
Om    = 2*pi*freqs ;

if printfreqs == 1
    digits(4) ;
    vpa(freqs)
end

if Nmodes > 8
    xax = linspace(0,Lx,Nx+1) ;
    yax = linspace(0,Ly,Ny+1) ;
    [X,Y] = meshgrid(xax,yax) ;
    colormap('parula') ;
    for m = 1 : 9
        mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
        subplot(3,3,m)
        mesh(X,Y,real(mdShape),(abs(mdShape)),'FaceColor','texturemap') ; view(2); axis equal; axis tight;
        xlabel('x (m)','interpreter','latex')
        ylabel('y (m)','interpreter','latex')
        tit = sprintf('Mode %d',m) ;
        title(tit) ;
    end
end

%------------ a bunch of functions ... these are only storing coefficients
%values really ....

    function [u00_c,u10_c,u20_c,u01_c,u02_c,u11_c] = D00_coeffs(K0y,R0y,Kx0,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4)


        u00_c = ((D3*((16*R0y*D3^3*hx^5 - 32*D1*nuy*D3^3*hx^4*nux + 32*D1*D3^3*hx^4 + 32*R0y*D3^3*hx^3*hy^2*nux - 32*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 32*D1*D3^3*hx^2*hy^2*nux + 24*R0y*D3^3*hx*hy^4*nux^2 - 16*D1*nuy*D3^3*hy^4*nux^3 + 16*D1*D3^3*hy^4*nux^2 + 8*R0y*D3^2*Rx0*hx^5*hy - 16*D1*nuy*D3^2*Rx0*hx^4*hy*nux + 16*D1*D3^2*Rx0*hx^4*hy + 16*R0y*D3^2*Rx0*hx^3*hy^3*nux - 16*D1*nuy*D3^2*Rx0*hx^2*hy^3*nux^2 + 16*D1*D3^2*Rx0*hx^2*hy^3*nux + 12*R0y*D3^2*Rx0*hx*hy^5*nux^2 + 8*D1*D3^2*Rx0*hy^5*nux^2 + 8*Kx0*R0y*D3^2*hx^5*hy^3 - 16*D1*Kx0*nuy*D3^2*hx^4*hy^3*nux + 16*D1*Kx0*D3^2*hx^4*hy^3 + 16*D4*R0y*D3^2*hx^3*hy^2 - 32*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 32*D1*D4*D3^2*hx^2*hy^2 + 24*D4*R0y*D3^2*hx*hy^4*nux - 16*D1*D4*nuy*D3^2*hy^4*nux^2 + 16*D1*D4*D3^2*hy^4*nux + 8*Kx0*R0y*D3*Rx0*hx^5*hy^4 - 8*D1*Kx0*nuy*D3*Rx0*hx^4*hy^4*nux + 16*D1*Kx0*D3*Rx0*hx^4*hy^4 + 8*D4*R0y*D3*Rx0*hx^3*hy^3 - 16*D1*D4*nuy*D3*Rx0*hx^2*hy^3*nux + 16*D1*D4*D3*Rx0*hx^2*hy^3 + 12*D4*R0y*D3*Rx0*hx*hy^5*nux + 8*D1*D4*D3*Rx0*hy^5*nux + 2*Kx0*R0y*Rx0^2*hx^5*hy^5 + 4*D1*Kx0*Rx0^2*hx^4*hy^5)/(D3*hx^2*(2*D3*hx^2 + Rx0*hx^2*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*(8*D1*D3*hx + 4*D3*R0y*hx^2 + 4*D3*R0y*hy^2*nux - 8*D1*D3*hx*nux*nuy))/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 6))/hy^4 + (D1*((24*Rx0*D1^3*hx^4*hy*nuy^2 - 16*D3*nux*D1^3*hx^4*nuy^3 + 16*D3*D1^3*hx^4*nuy^2 + 32*Rx0*D1^3*hx^2*hy^3*nuy - 32*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 32*D3*D1^3*hx^2*hy^2*nuy + 16*Rx0*D1^3*hy^5 - 32*D3*nux*D1^3*hy^4*nuy + 32*D3*D1^3*hy^4 + 12*Rx0*D1^2*R0y*hx^5*hy*nuy^2 + 8*D3*D1^2*R0y*hx^5*nuy^2 + 16*Rx0*D1^2*R0y*hx^3*hy^3*nuy - 16*D3*nux*D1^2*R0y*hx^3*hy^2*nuy^2 + 16*D3*D1^2*R0y*hx^3*hy^2*nuy + 8*Rx0*D1^2*R0y*hx*hy^5 - 16*D3*nux*D1^2*R0y*hx*hy^4*nuy + 16*D3*D1^2*R0y*hx*hy^4 + 24*D4*Rx0*D1^2*hx^4*hy*nuy - 16*D3*D4*nux*D1^2*hx^4*nuy^2 + 16*D3*D4*D1^2*hx^4*nuy + 8*K0y*Rx0*D1^2*hx^3*hy^5 - 16*D3*K0y*nux*D1^2*hx^3*hy^4*nuy + 16*D3*K0y*D1^2*hx^3*hy^4 + 16*D4*Rx0*D1^2*hx^2*hy^3 - 32*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 32*D3*D4*D1^2*hx^2*hy^2 + 12*D4*Rx0*D1*R0y*hx^5*hy*nuy + 8*D3*D4*D1*R0y*hx^5*nuy + 8*K0y*Rx0*D1*R0y*hx^4*hy^5 - 8*D3*K0y*nux*D1*R0y*hx^4*hy^4*nuy + 16*D3*K0y*D1*R0y*hx^4*hy^4 + 8*D4*Rx0*D1*R0y*hx^3*hy^3 - 16*D3*D4*nux*D1*R0y*hx^3*hy^2*nuy + 16*D3*D4*D1*R0y*hx^3*hy^2 + 2*K0y*Rx0*R0y^2*hx^5*hy^5 + 4*D3*K0y*R0y^2*hx^5*hy^4)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*(8*D1*D3*hy + 4*D1*Rx0*hy^2 + 4*D1*Rx0*hx^2*nuy - 8*D1*D3*hy*nux*nuy))/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 6))/hx^4 - ((D2 + D4)*((2*(8*D1*D3*hx + 4*D3*R0y*hx^2 + 4*D3*R0y*hy^2*nux - 8*D1*D3*hx*nux*nuy))/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (2*(8*D1*D3*hy + 4*D1*Rx0*hy^2 + 4*D1*Rx0*hx^2*nuy - 8*D1*D3*hy*nux*nuy))/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (2*D1*Rx0*nuy*hx^4*hy + 4*D1*D3*nuy*hx^4 + 2*D3*R0y*nux*hx*hy^4 + 4*D1*D3*nux*hy^4)/(hx*(2*D1*hy^2 + R0y*hx*hy^2)*(2*D3*hx + Rx0*hx*hy)) + (2*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) + (2*D3*hy^2*nux)/(hx*(2*D3*hx + Rx0*hx*hy)) - 4))/(hx^2*hy^2)) ;

        u10_c = (((D2 + D4)*((4*D3*hx^2 + 4*D3*nux*hy^2)/(hx*(2*D3*hx + Rx0*hx*hy)) - (2*(R0y*Rx0*hx*hy^2 - 2*D1*Rx0*hy^2 - 4*D1*D3*hy + 2*D3*R0y*hx*hy + 4*D1*D3*hy*nux*nuy))/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (4*D3*R0y*hx^3*hy^2 + 8*D1*D3*hx^2*hy^2 + 4*D3*R0y*nux*hx*hy^4 + 8*D1*D3*nux*hy^4)/(hx*(2*D1*hy^2 + R0y*hx*hy^2)*(2*D3*hx + Rx0*hx*hy)) + (8*D3*R0y*hy^2*nux)/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - 2))/(hx^2*hy^2) - (D3*((32*R0y*D3^3*hx^3*hy^2*nux - 32*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 32*D1*D3^3*hx^2*hy^2*nux + 32*R0y*D3^3*hx*hy^4*nux^2 - 32*D1*nuy*D3^3*hy^4*nux^3 + 32*D1*D3^3*hy^4*nux^2 + 16*R0y*Rx0*D3^2*hx^3*hy^3*nux + 16*D4*R0y*D3^2*hx^3*hy^2 + 16*D1*Rx0*D3^2*hx^2*hy^3*nux - 32*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 32*D1*D4*D3^2*hx^2*hy^2 + 16*R0y*Rx0*D3^2*hx*hy^5*nux^2 + 32*D4*R0y*D3^2*hx*hy^4*nux + 16*D1*Rx0*D3^2*hy^5*nux^2 - 32*D1*D4*nuy*D3^2*hy^4*nux^2 + 32*D1*D4*D3^2*hy^4*nux + 8*D4*R0y*Rx0*D3*hx^3*hy^3 + 16*D1*D4*Rx0*D3*hx^2*hy^3 + 16*D4*R0y*Rx0*D3*hx*hy^5*nux + 16*D1*D4*Rx0*D3*hy^5*nux)/(D3*hx^2*(2*D3*hx^2 + Rx0*hx^2*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (16*D3*R0y*hy^2*nux)/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))))/hy^4 - (D1*((4*(R0y*Rx0*hx*hy^2 - 2*D1*Rx0*hy^2 - 4*D1*D3*hy + 2*D3*R0y*hx*hy + 4*D1*D3*hy*nux*nuy))/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (16*Rx0*D1^3*hx^2*hy^3*nuy - 32*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 32*D3*D1^3*hx^2*hy^2*nuy + 16*Rx0*D1^3*hy^5 - 32*D3*nux*D1^3*hy^4*nuy + 32*D3*D1^3*hy^4 + 8*R0y*Rx0*D1^2*hx^3*hy^3*nuy - 16*D3*R0y*nux*D1^2*hx^3*hy^2*nuy^2 + 16*D3*R0y*D1^2*hx^3*hy^2*nuy + 16*D4*Rx0*D1^2*hx^2*hy^3 - 32*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 32*D3*D4*D1^2*hx^2*hy^2 + 8*R0y*Rx0*D1^2*hx*hy^5 - 16*D3*R0y*nux*D1^2*hx*hy^4*nuy + 16*D3*R0y*D1^2*hx*hy^4 + 8*D4*R0y*Rx0*D1*hx^3*hy^3 - 16*D3*D4*R0y*nux*D1*hx^3*hy^2*nuy + 16*D3*D4*R0y*D1*hx^3*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 4))/hx^4) ;
        
        u20_c = ((D1*((8*D1^3*D3*hy^4 + 4*D1^3*Rx0*hy^5 + 2*D1*D3*R0y^2*hx^2*hy^4 + D1*R0y^2*Rx0*hx^2*hy^5 + 8*D1^2*D3*R0y*hx*hy^4 + 4*D1^2*R0y*Rx0*hx*hy^5 - 8*D1^3*D3*hy^4*nux*nuy - 4*D1^2*D3*R0y*hx*hy^4*nux*nuy)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 1))/hx^4 - (((4*D1*D3*hy^4*nux + 2*D3*R0y*hx*hy^4*nux)/(hx*(2*D1*hy^2 + R0y*hx*hy^2)*(2*D3*hx + Rx0*hx*hy)) + (2*D3*hy^2*nux)/(hx*(2*D3*hx + Rx0*hx*hy)))*(D2 + D4))/(hx^2*hy^2) + (16*D1*D3^3*hy^4*nux^2 + 8*D1*D3^2*Rx0*hy^5*nux^2 - 16*D1*D3^3*hy^4*nux^3*nuy + 8*D3^3*R0y*hx*hy^4*nux^2 + 16*D1*D3^2*D4*hy^4*nux + 8*D3^2*D4*R0y*hx*hy^4*nux - 16*D1*D3^2*D4*hy^4*nux^2*nuy + 4*D3^2*R0y*Rx0*hx*hy^5*nux^2 + 8*D1*D3*D4*Rx0*hy^5*nux + 4*D3*D4*R0y*Rx0*hx*hy^5*nux)/(hx^2*hy^4*(2*D3*hx^2 + Rx0*hx^2*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))) ;

        u01_c = (((D2 + D4)*((4*D1*nuy*hx^2 + 4*D1*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - (2*(R0y*Rx0*hx^2*hy - 2*D3*R0y*hx^2 - 4*D1*D3*hx + 2*D1*Rx0*hx*hy + 4*D1*D3*hx*nux*nuy))/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (4*D1*Rx0*nuy*hx^4*hy + 8*D1*D3*nuy*hx^4 + 4*D1*Rx0*hx^2*hy^3 + 8*D1*D3*hx^2*hy^2)/(hx*(2*D1*hy^2 + R0y*hx*hy^2)*(2*D3*hx + Rx0*hx*hy)) + (8*D1*Rx0*hx^2*nuy)/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - 2))/(hx^2*hy^2) - (D3*((4*(R0y*Rx0*hx^2*hy - 2*D3*R0y*hx^2 - 4*D1*D3*hx + 2*D1*Rx0*hx*hy + 4*D1*D3*hx*nux*nuy))/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (16*R0y*D3^3*hx^5 - 32*D1*nuy*D3^3*hx^4*nux + 32*D1*D3^3*hx^4 + 16*R0y*D3^3*hx^3*hy^2*nux - 32*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 32*D1*D3^3*hx^2*hy^2*nux + 8*R0y*Rx0*D3^2*hx^5*hy - 16*D1*Rx0*nuy*D3^2*hx^4*hy*nux + 16*D1*Rx0*D3^2*hx^4*hy + 8*R0y*Rx0*D3^2*hx^3*hy^3*nux + 16*D4*R0y*D3^2*hx^3*hy^2 - 16*D1*Rx0*nuy*D3^2*hx^2*hy^3*nux^2 + 16*D1*Rx0*D3^2*hx^2*hy^3*nux - 32*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 32*D1*D4*D3^2*hx^2*hy^2 + 8*D4*R0y*Rx0*D3*hx^3*hy^3 - 16*D1*D4*Rx0*nuy*D3*hx^2*hy^3*nux + 16*D1*D4*Rx0*D3*hx^2*hy^3)/(D3*hx^2*(2*D3*hx^2 + Rx0*hx^2*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 4))/hy^4 - (D1*((32*Rx0*D1^3*hx^4*hy*nuy^2 - 32*D3*nux*D1^3*hx^4*nuy^3 + 32*D3*D1^3*hx^4*nuy^2 + 32*Rx0*D1^3*hx^2*hy^3*nuy - 32*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 32*D3*D1^3*hx^2*hy^2*nuy + 16*R0y*Rx0*D1^2*hx^5*hy*nuy^2 + 16*D3*R0y*D1^2*hx^5*nuy^2 + 32*D4*Rx0*D1^2*hx^4*hy*nuy - 32*D3*D4*nux*D1^2*hx^4*nuy^2 + 32*D3*D4*D1^2*hx^4*nuy + 16*R0y*Rx0*D1^2*hx^3*hy^3*nuy + 16*D3*R0y*D1^2*hx^3*hy^2*nuy + 16*D4*Rx0*D1^2*hx^2*hy^3 - 32*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 32*D3*D4*D1^2*hx^2*hy^2 + 16*D4*R0y*Rx0*D1*hx^5*hy*nuy + 16*D3*D4*R0y*D1*hx^5*nuy + 8*D4*R0y*Rx0*D1*hx^3*hy^3 + 16*D3*D4*R0y*D1*hx^3*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (16*D1*Rx0*hx^2*nuy)/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))))/hx^4) ;

        u02_c = ((D3*((8*D1*D3^3*hx^4 + 4*D3^3*R0y*hx^5 + 2*D1*D3*Rx0^2*hx^4*hy^2 + D3*R0y*Rx0^2*hx^5*hy^2 + 8*D1*D3^2*Rx0*hx^4*hy + 4*D3^2*R0y*Rx0*hx^5*hy - 8*D1*D3^3*hx^4*nux*nuy - 4*D1*D3^2*Rx0*hx^4*hy*nux*nuy)/(D3*hx^2*(2*D3*hx^2 + Rx0*hx^2*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 1))/hy^4 - (((4*D1*D3*hx^4*nuy + 2*D1*Rx0*hx^4*hy*nuy)/(hx*(2*D1*hy^2 + R0y*hx*hy^2)*(2*D3*hx + Rx0*hx*hy)) + (2*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2))*(D2 + D4))/(hx^2*hy^2) + (16*D1^3*D3*hx^4*nuy^2 + 8*D1^2*D3*R0y*hx^5*nuy^2 - 16*D1^3*D3*hx^4*nux*nuy^3 + 8*D1^3*Rx0*hx^4*hy*nuy^2 + 16*D1^2*D3*D4*hx^4*nuy + 8*D1^2*D4*Rx0*hx^4*hy*nuy - 16*D1^2*D3*D4*hx^4*nux*nuy^2 + 4*D1^2*R0y*Rx0*hx^5*hy*nuy^2 + 8*D1*D3*D4*R0y*hx^5*nuy + 4*D1*D4*R0y*Rx0*hx^5*hy*nuy)/(hx^4*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))) ;

        u11_c = ((16*Rx0*D1^3*hx^2*hy^3*nuy - 32*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 32*D3*D1^3*hx^2*hy^2*nuy + 8*R0y*Rx0*D1^2*hx^3*hy^3*nuy + 16*D3*R0y*D1^2*hx^3*hy^2*nuy + 16*D4*Rx0*D1^2*hx^2*hy^3 - 32*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 32*D3*D4*D1^2*hx^2*hy^2 + 8*D4*R0y*Rx0*D1*hx^3*hy^3 + 16*D3*D4*R0y*D1*hx^3*hy^2)/(hx^4*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - ((D2 + D4)*((2*D1*hy^2 - R0y*hx*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) + (2*D3*hx^2 - Rx0*hx^2*hy)/(hx*(2*D3*hx + Rx0*hx*hy)) + (- R0y*Rx0*hx^3*hy^3 + 2*D3*R0y*hx^3*hy^2 + 2*D1*Rx0*hx^2*hy^3 + 12*D1*D3*hx^2*hy^2)/(hx*(2*D1*hy^2 + R0y*hx*hy^2)*(2*D3*hx + Rx0*hx*hy)) - 1))/(hx^2*hy^2) + (16*R0y*D3^3*hx^3*hy^2*nux - 32*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 32*D1*D3^3*hx^2*hy^2*nux + 8*R0y*Rx0*D3^2*hx^3*hy^3*nux + 16*D4*R0y*D3^2*hx^3*hy^2 + 16*D1*Rx0*D3^2*hx^2*hy^3*nux - 32*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 32*D1*D4*D3^2*hx^2*hy^2 + 8*D4*R0y*Rx0*D3*hx^3*hy^3 + 16*D1*D4*Rx0*D3*hx^2*hy^3)/(hx^2*hy^4*(2*D3*hx^2 + Rx0*hx^2*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))) ;

    end

    function [u01_c,u11_c,u21_c,u00_c,u02_c,u03_c,u12_c,u10_c] = D01_coeffs(K0y,R0y,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4)

        u01_c  = ((D1*((28*Rx0*D1^3*hx^4*hy*nuy^2 - 40*D3*nux*D1^3*hx^4*nuy^3 + 40*D3*D1^3*hx^4*nuy^2 + 32*Rx0*D1^3*hx^2*hy^3*nuy - 64*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 64*D3*D1^3*hx^2*hy^2*nuy + 16*Rx0*D1^3*hy^5 - 32*D3*nux*D1^3*hy^4*nuy + 32*D3*D1^3*hy^4 + 14*Rx0*D1^2*R0y*hx^5*hy*nuy^2 + 20*D3*D1^2*R0y*hx^5*nuy^2 + 16*Rx0*D1^2*R0y*hx^3*hy^3*nuy + 32*D3*D1^2*R0y*hx^3*hy^2*nuy + 8*Rx0*D1^2*R0y*hx*hy^5 + 16*D3*D1^2*R0y*hx*hy^4 + 28*D4*Rx0*D1^2*hx^4*hy*nuy - 40*D3*D4*nux*D1^2*hx^4*nuy^2 + 40*D3*D4*D1^2*hx^4*nuy + 8*K0y*Rx0*D1^2*hx^3*hy^5 - 16*D3*K0y*nux*D1^2*hx^3*hy^4*nuy + 16*D3*K0y*D1^2*hx^3*hy^4 + 16*D4*Rx0*D1^2*hx^2*hy^3 - 32*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 32*D3*D4*D1^2*hx^2*hy^2 + 14*D4*Rx0*D1*R0y*hx^5*hy*nuy + 20*D3*D4*D1*R0y*hx^5*nuy + 8*K0y*Rx0*D1*R0y*hx^4*hy^5 - 8*D3*K0y*nux*D1*R0y*hx^4*hy^4*nuy + 16*D3*K0y*D1*R0y*hx^4*hy^4 + 8*D4*Rx0*D1*R0y*hx^3*hy^3 + 16*D3*D4*D1*R0y*hx^3*hy^2 + 2*K0y*Rx0*R0y^2*hx^5*hy^5 + 4*D3*K0y*R0y^2*hx^5*hy^4)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*(4*D1*nuy*hx^2 + 4*D1*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) + 6))/hx^4 + (D3*((R0y*Rx0*hx^2*hy - 2*D3*R0y*hx^2 - 4*D1*D3*hx + 2*D1*Rx0*hx*hy + 4*D1*D3*hx*nux*nuy)/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 6))/hy^4 - ((D2 + D4)*((2*(4*D1*nuy*hx^2 + 4*D1*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) + (2*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) + (4*D1*Rx0*hx^2*nuy)/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - 4))/(hx^2*hy^2)) ;

        u11_c  = ((((2*(2*D1*hy^2 - R0y*hx*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) - 2)*(D2 + D4))/(hx^2*hy^2) - (D1*((16*Rx0*D1^3*hx^2*hy^3*nuy - 32*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 32*D3*D1^3*hx^2*hy^2*nuy + 16*Rx0*D1^3*hy^5 - 32*D3*nux*D1^3*hy^4*nuy + 32*D3*D1^3*hy^4 + 8*R0y*Rx0*D1^2*hx^3*hy^3*nuy + 16*D3*R0y*D1^2*hx^3*hy^2*nuy + 16*D4*Rx0*D1^2*hx^2*hy^3 - 32*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 32*D3*D4*D1^2*hx^2*hy^2 + 8*R0y*Rx0*D1^2*hx*hy^5 + 16*D3*R0y*D1^2*hx*hy^4 + 8*D4*R0y*Rx0*D1*hx^3*hy^3 + 16*D3*D4*R0y*D1*hx^3*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*(2*D1*hy^2 - R0y*hx*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) + 4))/hx^4) ;

        u21_c  = ((D1*((8*D1^3*D3*hy^4 + 4*D1^3*Rx0*hy^5 + 2*D1*D3*R0y^2*hx^2*hy^4 + D1*R0y^2*Rx0*hx^2*hy^5 + 8*D1^2*D3*R0y*hx*hy^4 + 4*D1^2*R0y*Rx0*hx*hy^5 - 8*D1^3*D3*hy^4*nux*nuy - 4*D1^2*D3*R0y*hx*hy^4*nux*nuy)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 1))/hx^4) ;

        u00_c  = ((D1*((8*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) - (16*Rx0*D1^3*hx^4*hy*nuy^2 - 16*D3*nux*D1^3*hx^4*nuy^3 + 16*D3*D1^3*hx^4*nuy^2 + 16*Rx0*D1^3*hx^2*hy^3*nuy - 32*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 32*D3*D1^3*hx^2*hy^2*nuy + 8*R0y*Rx0*D1^2*hx^5*hy*nuy^2 + 8*D3*R0y*D1^2*hx^5*nuy^2 + 16*D4*Rx0*D1^2*hx^4*hy*nuy - 16*D3*D4*nux*D1^2*hx^4*nuy^2 + 16*D3*D4*D1^2*hx^4*nuy + 8*R0y*Rx0*D1^2*hx^3*hy^3*nuy - 8*D3*R0y*nux*D1^2*hx^3*hy^2*nuy^2 + 16*D3*R0y*D1^2*hx^3*hy^2*nuy + 8*D4*Rx0*D1^2*hx^2*hy^3 - 16*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 16*D3*D4*D1^2*hx^2*hy^2 + 8*D4*R0y*Rx0*D1*hx^5*hy*nuy + 8*D3*D4*R0y*D1*hx^5*nuy + 4*D4*R0y*Rx0*D1*hx^3*hy^3 - 8*D3*D4*R0y*nux*D1*hx^3*hy^2*nuy + 8*D3*D4*R0y*D1*hx^3*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))))/hx^4 + (D3*((8*D1*D3*hx + 4*D3*R0y*hx^2 + 4*D3*R0y*hy^2*nux - 8*D1*D3*hx*nux*nuy)/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - 4))/hy^4 + ((D2 + D4)*((8*D1*D3*hy + 4*D1*Rx0*hy^2 + 4*D1*Rx0*hx^2*nuy - 8*D1*D3*hy*nux*nuy)/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (4*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) - 2))/(hx^2*hy^2)) ;

        u02_c  = ((D1*((8*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) - (16*Rx0*D1^3*hx^4*hy*nuy^2 - 32*D3*nux*D1^3*hx^4*nuy^3 + 32*D3*D1^3*hx^4*nuy^2 + 16*Rx0*D1^3*hx^2*hy^3*nuy - 32*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 32*D3*D1^3*hx^2*hy^2*nuy + 8*R0y*Rx0*D1^2*hx^5*hy*nuy^2 + 16*D3*R0y*D1^2*hx^5*nuy^2 + 16*D4*Rx0*D1^2*hx^4*hy*nuy - 32*D3*D4*nux*D1^2*hx^4*nuy^2 + 32*D3*D4*D1^2*hx^4*nuy + 8*R0y*Rx0*D1^2*hx^3*hy^3*nuy + 16*D3*R0y*D1^2*hx^3*hy^2*nuy + 8*D4*Rx0*D1^2*hx^2*hy^3 - 16*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 16*D3*D4*D1^2*hx^2*hy^2 + 8*D4*R0y*Rx0*D1*hx^5*hy*nuy + 16*D3*D4*R0y*D1*hx^5*nuy + 4*D4*R0y*Rx0*D1*hx^3*hy^3 + 8*D3*D4*R0y*D1*hx^3*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))))/hx^4 - (4*D3)/hy^4 + ((D2 + D4)*((4*D1*nuy*hx^2 + 4*D1*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) + (4*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) - 2))/(hx^2*hy^2)) ;

        u03_c  = (D3/hy^4 + (8*D1^3*D3*hx^4*nuy^2 + 4*D1^2*D3*R0y*hx^5*nuy^2 - 8*D1^3*D3*hx^4*nux*nuy^3 + 4*D1^3*Rx0*hx^4*hy*nuy^2 + 8*D1^2*D3*D4*hx^4*nuy + 4*D1^2*D4*Rx0*hx^4*hy*nuy - 8*D1^2*D3*D4*hx^4*nux*nuy^2 + 2*D1^2*R0y*Rx0*hx^5*hy*nuy^2 + 4*D1*D3*D4*R0y*hx^5*nuy + 2*D1*D4*R0y*Rx0*hx^5*hy*nuy)/(hx^4*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (2*D1*nuy*(D2 + D4))/(hy^2*(2*D1*hy^2 + R0y*hx*hy^2))) ;

        u12_c  = ((8*Rx0*D1^3*hx^2*hy^3*nuy - 16*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 16*D3*D1^3*hx^2*hy^2*nuy + 4*R0y*Rx0*D1^2*hx^3*hy^3*nuy + 8*D3*R0y*D1^2*hx^3*hy^2*nuy + 8*D4*Rx0*D1^2*hx^2*hy^3 - 16*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 16*D3*D4*D1^2*hx^2*hy^2 + 4*D4*R0y*Rx0*D1*hx^3*hy^3 + 8*D3*D4*R0y*D1*hx^3*hy^2)/(hx^4*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (((2*D1*hy^2 - R0y*hx*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - 1)*(D2 + D4))/(hx^2*hy^2)) ;

        u10_c  = ((((R0y*Rx0*hx*hy^2 - 2*D1*Rx0*hy^2 - 4*D1*D3*hy + 2*D3*R0y*hx*hy + 4*D1*D3*hy*nux*nuy)/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 1)*(D2 + D4))/(hx^2*hy^2) + (8*Rx0*D1^3*hx^2*hy^3*nuy - 16*D3*nux*D1^3*hx^2*hy^2*nuy^2 + 16*D3*D1^3*hx^2*hy^2*nuy + 4*R0y*Rx0*D1^2*hx^3*hy^3*nuy - 8*D3*R0y*nux*D1^2*hx^3*hy^2*nuy^2 + 8*D3*R0y*D1^2*hx^3*hy^2*nuy + 8*D4*Rx0*D1^2*hx^2*hy^3 - 16*D3*D4*nux*D1^2*hx^2*hy^2*nuy + 16*D3*D4*D1^2*hx^2*hy^2 + 4*D4*R0y*Rx0*D1*hx^3*hy^3 - 8*D3*D4*R0y*nux*D1*hx^3*hy^2*nuy + 8*D3*D4*R0y*D1*hx^3*hy^2)/(hx^4*hy^4*(2*D1 + R0y*hx)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*D3^2*R0y*nux)/(hx*hy^2*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))) ;

    end

    function [u02_c,u12_c,u22_c,u01_c,u03_c,u04_c,u00_c,u13_c,u11_c] = D02_coeffs(K0y,R0y,hx,hy,nuy,D1,D2,D3,D4)

        u02_c = ((6*D3)/hy^4 + (D1*((12*D1^2*hx^4*nuy^2 + 16*D1^2*hx^2*hy^2*nuy + 8*D1^2*hy^4 + 12*D4*D1*hx^4*nuy + 4*K0y*D1*hx^3*hy^4 + 8*D4*D1*hx^2*hy^2 + 2*K0y*R0y*hx^4*hy^4)/(D1*hy^4*(2*D1 + R0y*hx)) - (4*(4*D1*nuy*hx^2 + 4*D1*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) + 6))/hx^4 - ((D2 + D4)*((2*(4*D1*nuy*hx^2 + 4*D1*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) + (4*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) - 4))/(hx^2*hy^2)) ;

        u12_c = ((((2*(2*D1*hy^2 - R0y*hx*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) - 2)*(D2 + D4))/(hx^2*hy^2) - (D1*((8*nuy*D1^2*hx^2*hy^2 + 8*D1^2*hy^4 + 8*D4*D1*hx^2*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)) - (4*(2*D1*hy^2 - R0y*hx*hy^2))/(2*D1*hy^2 + R0y*hx*hy^2) + 4))/hx^4) ;

        u22_c = ((D1*((2*D1^2*hy^4 + R0y*hx*D1*hy^4)/(D1*hy^4*(2*D1 + R0y*hx)) + 1))/hx^4) ;

        u01_c = (((D2 + D4)*((4*D1*nuy*hx^2 + 4*D1*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) + (4*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) - 2))/(hx^2*hy^2) - (D1*((8*D1^2*hx^4*nuy^2 + 8*D1^2*hx^2*hy^2*nuy + 8*D4*D1*hx^4*nuy + 4*D4*D1*hx^2*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)) - (8*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2)))/hx^4 - (4*D3)/hy^4) ;

        u03_c = (((D2 + D4)*((4*D1*nuy*hx^2 + 4*D1*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) + (4*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2) - 2))/(hx^2*hy^2) - (D1*((8*D1^2*hx^4*nuy^2 + 8*D1^2*hx^2*hy^2*nuy + 8*D4*D1*hx^4*nuy + 4*D4*D1*hx^2*hy^2)/(D1*hy^4*(2*D1 + R0y*hx)) - (8*D1*hx^2*nuy)/(2*D1*hy^2 + R0y*hx*hy^2)))/hx^4 - (4*D3)/hy^4) ;

        u04_c = (D3/hy^4 + (2*D1^2*hx^4*nuy^2 + 2*D4*D1*hx^4*nuy)/(hx^4*hy^4*(2*D1 + R0y*hx)) - (2*D1*nuy*(D2 + D4))/(hy^2*(2*D1*hy^2 + R0y*hx*hy^2))) ;

        u00_c = (D3/hy^4 + (2*D1^2*hx^4*nuy^2 + 2*D4*D1*hx^4*nuy)/(hx^4*hy^4*(2*D1 + R0y*hx)) - (2*D1*nuy*(D2 + D4))/(hy^2*(2*D1*hy^2 + R0y*hx*hy^2))) ;

        u13_c = ((4*nuy*D1^2*hx^2*hy^2 + 4*D4*D1*hx^2*hy^2)/(hx^4*hy^4*(2*D1 + R0y*hx)) - (((2*D1*hy^2 - R0y*hx*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - 1)*(D2 + D4))/(hx^2*hy^2)) ;

        u11_c = ((4*nuy*D1^2*hx^2*hy^2 + 4*D4*D1*hx^2*hy^2)/(hx^4*hy^4*(2*D1 + R0y*hx)) - (((2*D1*hy^2 - R0y*hx*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - 1)*(D2 + D4))/(hx^2*hy^2)) ;

    end

    function [u10_c, u20_c, u30_c, u00_c, u11_c, u12_c, u21_c, u01_c] = D10_coeffs(R0y,Kx0,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4)

        u10_c = ((D3*((16*R0y*D3^3*hx^5 - 32*D1*nuy*D3^3*hx^4*nux + 32*D1*D3^3*hx^4 + 32*R0y*D3^3*hx^3*hy^2*nux - 64*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 64*D1*D3^3*hx^2*hy^2*nux + 28*R0y*D3^3*hx*hy^4*nux^2 - 40*D1*nuy*D3^3*hy^4*nux^3 + 40*D1*D3^3*hy^4*nux^2 + 8*R0y*D3^2*Rx0*hx^5*hy + 16*D1*D3^2*Rx0*hx^4*hy + 16*R0y*D3^2*Rx0*hx^3*hy^3*nux + 32*D1*D3^2*Rx0*hx^2*hy^3*nux + 14*R0y*D3^2*Rx0*hx*hy^5*nux^2 + 20*D1*D3^2*Rx0*hy^5*nux^2 + 8*Kx0*R0y*D3^2*hx^5*hy^3 - 16*D1*Kx0*nuy*D3^2*hx^4*hy^3*nux + 16*D1*Kx0*D3^2*hx^4*hy^3 + 16*D4*R0y*D3^2*hx^3*hy^2 - 32*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 32*D1*D4*D3^2*hx^2*hy^2 + 28*D4*R0y*D3^2*hx*hy^4*nux - 40*D1*D4*nuy*D3^2*hy^4*nux^2 + 40*D1*D4*D3^2*hy^4*nux + 8*Kx0*R0y*D3*Rx0*hx^5*hy^4 - 8*D1*Kx0*nuy*D3*Rx0*hx^4*hy^4*nux + 16*D1*Kx0*D3*Rx0*hx^4*hy^4 + 8*D4*R0y*D3*Rx0*hx^3*hy^3 + 16*D1*D4*D3*Rx0*hx^2*hy^3 + 14*D4*R0y*D3*Rx0*hx*hy^5*nux + 20*D1*D4*D3*Rx0*hy^5*nux + 2*Kx0*R0y*Rx0^2*hx^5*hy^5 + 4*D1*Kx0*Rx0^2*hx^4*hy^5)/(D3*hx^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*(4*D3*hx^2 + 4*D3*nux*hy^2))/(hx*(2*D3*hx + Rx0*hx*hy)) + 6))/hy^4 + (D1*((R0y*Rx0*hx*hy^2 - 2*D1*Rx0*hy^2 - 4*D1*D3*hy + 2*D3*R0y*hx*hy + 4*D1*D3*hy*nux*nuy)/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 6))/hx^4 - ((D2 + D4)*((2*(4*D3*hx^2 + 4*D3*nux*hy^2))/(hx*(2*D3*hx + Rx0*hx*hy)) + (2*D3*hy^2*nux)/(2*D3*hx^2 + Rx0*hx^2*hy) + (4*D3*R0y*hy^2*nux)/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - 4))/(hx^2*hy^2)) ;

        u20_c = (((D2 + D4)*((4*D3*hx^2 + 4*D3*nux*hy^2)/(2*D3*hx^2 + Rx0*hx^2*hy) + (4*D3*hy^2*nux)/(hx*(2*D3*hx + Rx0*hx*hy)) - 2))/(hx^2*hy^2) - (D3*((16*R0y*D3^3*hx^3*hy^2*nux - 32*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 32*D1*D3^3*hx^2*hy^2*nux + 16*R0y*D3^3*hx*hy^4*nux^2 - 32*D1*nuy*D3^3*hy^4*nux^3 + 32*D1*D3^3*hy^4*nux^2 + 8*R0y*Rx0*D3^2*hx^3*hy^3*nux + 8*D4*R0y*D3^2*hx^3*hy^2 + 16*D1*Rx0*D3^2*hx^2*hy^3*nux - 16*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 16*D1*D4*D3^2*hx^2*hy^2 + 8*R0y*Rx0*D3^2*hx*hy^5*nux^2 + 16*D4*R0y*D3^2*hx*hy^4*nux + 16*D1*Rx0*D3^2*hy^5*nux^2 - 32*D1*D4*nuy*D3^2*hy^4*nux^2 + 32*D1*D4*D3^2*hy^4*nux + 4*D4*R0y*Rx0*D3*hx^3*hy^3 + 8*D1*D4*Rx0*D3*hx^2*hy^3 + 8*D4*R0y*Rx0*D3*hx*hy^5*nux + 16*D1*D4*Rx0*D3*hy^5*nux)/(D3*hx^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (8*D3*hy^2*nux)/(hx*(2*D3*hx + Rx0*hx*hy))))/hy^4 - (4*D1)/hx^4) ;

        u30_c = (D1/hx^4 + (8*D1*D3^3*hy^4*nux^2 + 4*D1*D3^2*Rx0*hy^5*nux^2 - 8*D1*D3^3*hy^4*nux^3*nuy + 4*D3^3*R0y*hx*hy^4*nux^2 + 8*D1*D3^2*D4*hy^4*nux + 4*D3^2*D4*R0y*hx*hy^4*nux - 8*D1*D3^2*D4*hy^4*nux^2*nuy + 2*D3^2*R0y*Rx0*hx*hy^5*nux^2 + 4*D1*D3*D4*Rx0*hy^5*nux + 2*D3*D4*R0y*Rx0*hx*hy^5*nux)/(hx^4*hy^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (2*D3*nux*(D2 + D4))/(hx^2*(2*D3*hx^2 + Rx0*hx^2*hy))) ;

        u00_c = ((D3*((8*D3*hy^2*nux)/(hx*(2*D3*hx + Rx0*hx*hy)) - (16*R0y*D3^3*hx^3*hy^2*nux - 32*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 32*D1*D3^3*hx^2*hy^2*nux + 16*R0y*D3^3*hx*hy^4*nux^2 - 16*D1*nuy*D3^3*hy^4*nux^3 + 16*D1*D3^3*hy^4*nux^2 + 8*R0y*Rx0*D3^2*hx^3*hy^3*nux + 8*D4*R0y*D3^2*hx^3*hy^2 - 8*D1*Rx0*nuy*D3^2*hx^2*hy^3*nux^2 + 16*D1*Rx0*D3^2*hx^2*hy^3*nux - 16*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 16*D1*D4*D3^2*hx^2*hy^2 + 8*R0y*Rx0*D3^2*hx*hy^5*nux^2 + 16*D4*R0y*D3^2*hx*hy^4*nux + 8*D1*Rx0*D3^2*hy^5*nux^2 - 16*D1*D4*nuy*D3^2*hy^4*nux^2 + 16*D1*D4*D3^2*hy^4*nux + 4*D4*R0y*Rx0*D3*hx^3*hy^3 - 8*D1*D4*Rx0*nuy*D3*hx^2*hy^3*nux + 8*D1*D4*Rx0*D3*hx^2*hy^3 + 8*D4*R0y*Rx0*D3*hx*hy^5*nux + 8*D1*D4*Rx0*D3*hy^5*nux)/(D3*hx^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))))/hy^4 + (D1*((8*D1*D3*hy + 4*D1*Rx0*hy^2 + 4*D1*Rx0*hx^2*nuy - 8*D1*D3*hy*nux*nuy)/(hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - 4))/hx^4 + ((D2 + D4)*((8*D1*D3*hx + 4*D3*R0y*hx^2 + 4*D3*R0y*hy^2*nux - 8*D1*D3*hx*nux*nuy)/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + (4*D3*hy^2*nux)/(hx*(2*D3*hx + Rx0*hx*hy)) - 2))/(hx^2*hy^2)) ;

        u11_c = ((((2*(2*D3*hx^2 - Rx0*hx^2*hy))/(hx*(2*D3*hx + Rx0*hx*hy)) - 2)*(D2 + D4))/(hx^2*hy^2) - (D3*((16*R0y*D3^3*hx^5 - 32*D1*nuy*D3^3*hx^4*nux + 32*D1*D3^3*hx^4 + 16*R0y*D3^3*hx^3*hy^2*nux - 32*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 32*D1*D3^3*hx^2*hy^2*nux + 8*R0y*Rx0*D3^2*hx^5*hy + 16*D1*Rx0*D3^2*hx^4*hy + 8*R0y*Rx0*D3^2*hx^3*hy^3*nux + 16*D4*R0y*D3^2*hx^3*hy^2 + 16*D1*Rx0*D3^2*hx^2*hy^3*nux - 32*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 32*D1*D4*D3^2*hx^2*hy^2 + 8*D4*R0y*Rx0*D3*hx^3*hy^3 + 16*D1*D4*Rx0*D3*hx^2*hy^3)/(D3*hx^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*(2*D3*hx^2 - Rx0*hx^2*hy))/(hx*(2*D3*hx + Rx0*hx*hy)) + 4))/hy^4) ;

        u12_c = ((D3*((8*D1*D3^3*hx^4 + 4*D3^3*R0y*hx^5 + 2*D1*D3*Rx0^2*hx^4*hy^2 + D3*R0y*Rx0^2*hx^5*hy^2 + 8*D1*D3^2*Rx0*hx^4*hy + 4*D3^2*R0y*Rx0*hx^5*hy - 8*D1*D3^3*hx^4*nux*nuy - 4*D1*D3^2*Rx0*hx^4*hy*nux*nuy)/(D3*hx^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 1))/hy^4) ;

        u21_c = ((8*R0y*D3^3*hx^3*hy^2*nux - 16*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 16*D1*D3^3*hx^2*hy^2*nux + 4*R0y*Rx0*D3^2*hx^3*hy^3*nux + 8*D4*R0y*D3^2*hx^3*hy^2 + 8*D1*Rx0*D3^2*hx^2*hy^3*nux - 16*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 16*D1*D4*D3^2*hx^2*hy^2 + 4*D4*R0y*Rx0*D3*hx^3*hy^3 + 8*D1*D4*Rx0*D3*hx^2*hy^3)/(hx^4*hy^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (((2*D3*hx^2 - Rx0*hx^2*hy)/(2*D3*hx^2 + Rx0*hx^2*hy) - 1)*(D2 + D4))/(hx^2*hy^2)) ;

        u01_c = ((((R0y*Rx0*hx^2*hy - 2*D3*R0y*hx^2 - 4*D1*D3*hx + 2*D1*Rx0*hx*hy + 4*D1*D3*hx*nux*nuy)/(hx*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) + 1)*(D2 + D4))/(hx^2*hy^2) + (8*R0y*D3^3*hx^3*hy^2*nux - 16*D1*nuy*D3^3*hx^2*hy^2*nux^2 + 16*D1*D3^3*hx^2*hy^2*nux + 4*R0y*Rx0*D3^2*hx^3*hy^3*nux + 8*D4*R0y*D3^2*hx^3*hy^2 - 8*D1*Rx0*nuy*D3^2*hx^2*hy^3*nux^2 + 8*D1*Rx0*D3^2*hx^2*hy^3*nux - 16*D1*D4*nuy*D3^2*hx^2*hy^2*nux + 16*D1*D4*D3^2*hx^2*hy^2 + 4*D4*R0y*Rx0*D3*hx^3*hy^3 - 8*D1*D4*Rx0*nuy*D3*hx^2*hy^3*nux + 8*D1*D4*Rx0*D3*hx^2*hy^3)/(hx^4*hy^4*(2*D3 + Rx0*hy)*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy)) - (4*D1^2*Rx0*nuy)/(hx^2*hy*(4*D1*D3 + 2*D3*R0y*hx + 2*D1*Rx0*hy - 4*D1*D3*nux*nuy + R0y*Rx0*hx*hy))) ;

    end

    function [u11_c,u12_c,u13_c,u10_c,u01_c,u21_c,u31_c,u22_c,u20_c,u00_c,u02_c] = D11_coeffs(R0y,Rx0,hx,hy,nux,nuy,D1,D2,D3,D4)

        u11_c = ((4*(D2 + D4))/(hx^2*hy^2) - (D1*((2*D1*hy^2 - R0y*hx*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - 6))/hx^4 - (D3*((2*D3*hx^2 - Rx0*hx^2*hy)/(hx*(2*D3*hx + Rx0*hx*hy)) - 6))/hy^4) ;

        u12_c = (- (4*D3)/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u13_c =  (D3/hy^4) ;

        u10_c = ((D3*((4*D3*hx^2 + 4*D3*nux*hy^2)/(hx*(2*D3*hx + Rx0*hx*hy)) - 4))/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u01_c = ((D1*((4*D1*nuy*hx^2 + 4*D1*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - 4))/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u21_c = (- (4*D1)/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u31_c = (D1/hx^4) ;

        u22_c = ((D2 + D4)/(hx^2*hy^2)) ;

        u20_c = ((D2 + D4)/(hx^2*hy^2) - (2*D3^2*nux)/(hx*hy^2*(2*D3*hx + Rx0*hx*hy))) ;

        u00_c = ((D2 + D4)/(hx^2*hy^2) - (2*D1^2*nuy)/(hx^2*(2*D1*hy^2 + R0y*hx*hy^2)) - (2*D3^2*nux)/(hx*hy^2*(2*D3*hx + Rx0*hx*hy))) ;

        u02_c = ((D2 + D4)/(hx^2*hy^2) - (2*D1^2*nuy)/(hx^2*(2*D1*hy^2 + R0y*hx*hy^2))) ;

    end

    function [u12_c,u13_c,u14_c,u11_c,u10_c,u02_c,u22_c,u32_c,u23_c,u21_c,u01_c,u03_c] = D12_coeffs(R0y,hx,hy,nuy,D1,D2,D3,D4)

        u12_c = ((6*D3)/hy^4 + (4*(D2 + D4))/(hx^2*hy^2) - (D1*((2*D1*hy^2 - R0y*hx*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - 6))/hx^4) ;

        u13_c = (- (4*D3)/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u14_c = (D3/hy^4) ;

        u11_c = (- (4*D3)/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u10_c = (D3/hy^4) ;

        u02_c = ((D1*((4*D1*nuy*hx^2 + 4*D1*hy^2)/(2*D1*hy^2 + R0y*hx*hy^2) - 4))/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u22_c = (- (4*D1)/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u32_c = (D1/hx^4) ;

        u23_c = ((D2 + D4)/(hx^2*hy^2)) ;

        u21_c = ((D2 + D4)/(hx^2*hy^2)) ;

        u01_c = ((D2 + D4)/(hx^2*hy^2) - (2*D1^2*nuy)/(hx^2*(2*D1*hy^2 + R0y*hx*hy^2))) ;

        u03_c = ((D2 + D4)/(hx^2*hy^2) - (2*D1^2*nuy)/(hx^2*(2*D1*hy^2 + R0y*hx*hy^2))) ;

    end

    function [u20_c,u21_c,u22_c,u10_c,u30_c,u40_c,u00_c,u31_c,u11_c] = D20_coeffs(Kx0,Rx0,hx,hy,nux,D1,D2,D3,D4)

        u20_c = ((6*D1)/hx^4 + (D3*((8*D3^2*hx^4 + 16*D3^2*hx^2*hy^2*nux + 12*D3^2*hy^4*nux^2 + 4*Kx0*D3*hx^4*hy^3 + 8*D4*D3*hx^2*hy^2 + 12*D4*D3*hy^4*nux + 2*Kx0*Rx0*hx^4*hy^4)/(D3*hx^4*(2*D3 + Rx0*hy)) - (4*(4*D3*hx^2 + 4*D3*nux*hy^2))/(2*D3*hx^2 + Rx0*hx^2*hy) + 6))/hy^4 - ((D2 + D4)*((2*(4*D3*hx^2 + 4*D3*nux*hy^2))/(2*D3*hx^2 + Rx0*hx^2*hy) + (2*D3*hy^2*nux)/(2*D3*hx^2 + Rx0*hx^2*hy) + (2*D3*hy^2*nux)/(hx*(2*D3*hx + Rx0*hx*hy)) - 4))/(hx^2*hy^2)) ; 

        u21_c = ((((2*(2*D3*hx^2 - Rx0*hx^2*hy))/(2*D3*hx^2 + Rx0*hx^2*hy) - 2)*(D2 + D4))/(hx^2*hy^2) - (D3*((8*D3^2*hx^4 + 8*nux*D3^2*hx^2*hy^2 + 8*D4*D3*hx^2*hy^2)/(D3*hx^4*(2*D3 + Rx0*hy)) - (4*(2*D3*hx^2 - Rx0*hx^2*hy))/(2*D3*hx^2 + Rx0*hx^2*hy) + 4))/hy^4) ;

        u22_c = ((D3*((2*D3^2*hx^4 + Rx0*hy*D3*hx^4)/(D3*hx^4*(2*D3 + Rx0*hy)) + 1))/hy^4) ;

        u10_c = (((D2 + D4)*((4*D3*hx^2 + 4*D3*nux*hy^2)/(hx*(2*D3*hx + Rx0*hx*hy)) + (4*D3*hy^2*nux)/(2*D3*hx^2 + Rx0*hx^2*hy) - 2))/(hx^2*hy^2) - (D3*((8*D3^2*hx^2*hy^2*nux + 8*D3^2*hy^4*nux^2 + 4*D4*D3*hx^2*hy^2 + 8*D4*D3*hy^4*nux)/(D3*hx^4*(2*D3 + Rx0*hy)) - (8*D3*hy^2*nux)/(2*D3*hx^2 + Rx0*hx^2*hy)))/hy^4 - (4*D1)/hx^4) ;

        u30_c = (((D2 + D4)*((4*D3*hx^2 + 4*D3*nux*hy^2)/(2*D3*hx^2 + Rx0*hx^2*hy) + (4*D3*hy^2*nux)/(2*D3*hx^2 + Rx0*hx^2*hy) - 2))/(hx^2*hy^2) - (D3*((8*D3^2*hx^2*hy^2*nux + 8*D3^2*hy^4*nux^2 + 4*D4*D3*hx^2*hy^2 + 8*D4*D3*hy^4*nux)/(D3*hx^4*(2*D3 + Rx0*hy)) - (8*D3*hy^2*nux)/(2*D3*hx^2 + Rx0*hx^2*hy)))/hy^4 - (4*D1)/hx^4) ; 

        u40_c = (D1/hx^4 + (2*D3^2*hy^4*nux^2 + 2*D4*D3*hy^4*nux)/(hx^4*hy^4*(2*D3 + Rx0*hy)) - (2*D3*nux*(D2 + D4))/(hx^2*(2*D3*hx^2 + Rx0*hx^2*hy))) ;

        u00_c = (D1/hx^4 + (2*D3^2*hy^4*nux^2 + 2*D4*D3*hy^4*nux)/(hx^4*hy^4*(2*D3 + Rx0*hy)) - (2*D3*nux*(D2 + D4))/(hx^3*(2*D3*hx + Rx0*hx*hy))) ;

        u31_c = ((4*nux*D3^2*hx^2*hy^2 + 4*D4*D3*hx^2*hy^2)/(hx^4*hy^4*(2*D3 + Rx0*hy)) - (((2*D3*hx^2 - Rx0*hx^2*hy)/(2*D3*hx^2 + Rx0*hx^2*hy) - 1)*(D2 + D4))/(hx^2*hy^2)) ;

        u11_c = ((4*nux*D3^2*hx^2*hy^2 + 4*D4*D3*hx^2*hy^2)/(hx^4*hy^4*(2*D3 + Rx0*hy)) - (((2*D3*hx^2 - Rx0*hx^2*hy)/(hx*(2*D3*hx + Rx0*hx*hy)) - 1)*(D2 + D4))/(hx^2*hy^2)) ; 

    end

    function [u21_c,u22_c,u23_c,u20_c,u11_c,u31_c,u41_c,u01_c,u32_c,u30_c,u10_c,u12_c] = D21_coeffs(Rx0,hx,hy,nux,D1,D2,D3,D4)

        u21_c = ((6*D1)/hx^4 + (4*(D2 + D4))/(hx^2*hy^2) - (D3*((2*D3*hx^2 - Rx0*hx^2*hy)/(2*D3*hx^2 + Rx0*hx^2*hy) - 6))/hy^4) ; 

        u22_c = (- (4*D3)/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u23_c = (D3/hy^4) ;

        u20_c = ((D3*((4*D3*hx^2 + 4*D3*nux*hy^2)/(2*D3*hx^2 + Rx0*hx^2*hy) - 4))/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u11_c = (- (4*D1)/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u31_c = (- (4*D1)/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u41_c = (D1/hx^4) ; 

        u01_c = (D1/hx^4) ;

        u32_c = ((D2 + D4)/(hx^2*hy^2)) ;

        u30_c = ((D2 + D4)/(hx^2*hy^2) - (2*D3^2*nux)/(hy^2*(2*D3*hx^2 + Rx0*hx^2*hy))) ;

        u10_c = ((D2 + D4)/(hx^2*hy^2) - (2*D3^2*nux)/(hy^2*(2*D3*hx^2 + Rx0*hx^2*hy))) ;

        u12_c = ((D2 + D4)/(hx^2*hy^2)) ;

    end


    function [u20_c,u11_c,u21_c,u31_c,u02_c,u12_c,u22_c,u32_c,u42_c,u13_c,u23_c,u33_c,u24_c] = D22_coeffs(hx,hy,D1,D2,D3,D4)

        u20_c = (D3)/hy^4 ;

        u11_c = ((D2 + D4)/(hx^2*hy^2)) ;

        u21_c = (- (4*D3)/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u31_c = ((D2 + D4)/(hx^2*hy^2)) ;

        u02_c = (D1/hx^4) ;

        u12_c = (- (4*D1)/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u22_c = ((6*D1)/hx^4 + (6*D3)/hy^4 + (4*(D2 + D4))/(hx^2*hy^2)) ; 

        u32_c = (- (4*D1)/hx^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u42_c = (D1/hx^4) ;

        u13_c = ((D2 + D4)/(hx^2*hy^2)) ;

        u23_c = (- (4*D3)/hy^4 - (2*(D2 + D4))/(hx^2*hy^2)) ;

        u33_c = ((D2 + D4)/(hx^2*hy^2)) ; 

        u24_c = (D3/hy^4) ;

    end



end
