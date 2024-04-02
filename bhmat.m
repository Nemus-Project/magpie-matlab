function biHarm = bhmat(BCs,Nxy,h,Lz,E,nu)
% biHarm = BHMAT(BCs,Nxy,h,D,nu) compute a biharmonic for a Thin Plate
% of size [Nx Ny]
%
%   Arguments:
%       BCS          %-- density [kg/m^3]
%       Nxy          %-- Number of grid points [Nx Ny]
%       h            %-- Grid spacing
%       Lz           %-- plate thickness
%       E            %-- Young's mod [Pa]
%       nu           %-- poisson's ratio
%
%% validate
validateattributes(BCs,      {'double'}, {'size', [4,2]});
validateattributes(Nxy,      {'double'}, {'numel', 2, 'positive','integer'});
% validateattributes(h,        {'double'}, {'nonempty'});
% validateattributes(Lz,       {'double'}, {'nonempty'});
% validateattributes(E,        {'double'}, {'nonempty'});
% validateattributes(nu,       {'double'}, {'nonempty'});

%% Unpack Variables
pack_BCs = num2cell(BCs);
[K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = pack_BCs{:};
Nx = Nxy(1);
Ny = Nxy(2);

D = E * Lz^3 / 12 / (1-nu^2);
%% MATRIX BUILDER

%%--- build matrix in blocks

a0 = ones(Ny+1,1) ;
a1 = ones(Ny,1) ;
a2 = ones(Ny-1,1) ;

[D00u00,D00u10,D00u20,D00u01,D00u02,D00u11] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[D01u01,D01u11,D01u21,D01u00,D01u02,D01u03,D01u12,D01u10] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[D02u02,D02u12,D02u22,D02u01,D02u03,D02u04,D02u00,D02u13,D02u11] = D02_coeffs(K0y,R0y,h,D,nu) ;
[D0Nu0N,D0Nu1N,D0Nu2N,D0Nu0Nm1,D0Nu0Nm2,D0Nu1Nm1] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,D0Nm1u1Nm1,D0Nm1u2Nm1,D0Nm1u0N,D0Nm1u0Nm2,D0Nm1u0Nm3,D0Nm1u1Nm2,D0Nm1u1N] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;
%%
%%-- Blk11
D0 = D02u02*a0; D1 = D02u03*a1; D2 = D02u04*a2 ; Dm1 = D02u01*a1; Dm2 = D02u00*a2 ;

Blk11               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
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

%%
%%%-- Blk12
D0 = D02u12*a0 ; D1 = D02u13*a1 ; Dm1 = D02u11*a1 ;

Blk12               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
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

%%
%%-- Blk13
D0 = D02u22*a0 ;

Blk13               = sparse(diag(D0))   ;
Blk13(1,1)          = D00u20 ;
Blk13(2,2)          = D01u21 ;
Blk13(end,end)      = D0Nu2N ;
Blk13(end-1,end-1)  = D0Nm1u2Nm1 ;

[D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[D11u11,D11u12,D11u13,D11u10,D11u01,D11u21,D11u31,D11u22,D11u20,D11u00,D11u02] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[D12u12,D12u13,D12u14,D12u11,D12u10,D12u02,D12u22,D12u32,D12u23,D12u21,D12u01,D12u03] = D12_coeffs(R0y,h,D,nu) ;
[D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[D1Nm1u1Nm1,D1Nm1u1Nm2,D1Nm1u1Nm3,D1Nm1u1N,D1Nm1u0Nm1,D1Nm1u2Nm1,D1Nm1u3Nm1,D1Nm1u2Nm2,D1Nm1u2N,D1Nm1u0N,D1Nm1u0Nm2] = D11_coeffs(R0y,RxL,h,D,nu) ;
%%
%%-- Blk21
D0 = D12u02*a0 ; D1 = D12u03*a1 ; Dm1 = D12u01*a1 ;

Blk21               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
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

%%
%%-- Blk22
D0 = D12u12*a0 ; D1 = D12u13*a1 ; D2 = D12u14*a2 ; Dm1 = D12u11*a1 ; Dm2 = D12u10*a2 ;

Blk22               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
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

%%
%%-- Blk23
D0 = D12u22*a0 ; D1 = D12u23*a1 ; Dm1 = D12u21*a1 ;

Blk23               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
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

%%
%%-- Blk24
D0 = D12u32*a0 ;

Blk24               = sparse(diag(D0))   ;
Blk24(1,1)          = D10u30 ;
Blk24(2,2)          = D11u31 ;
Blk24(end,end)      = D1Nu3N ;
Blk24(end-1,end-1)  = D1Nm1u3Nm1 ;

[D20u20,D20u21,D20u22,D20u10,D20u30,D20u40,D20u00,D20u31,D20u11] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[D21u21,D21u22,D21u23,D21u20,D21u11,D21u31,D21u41,D21u01,D21u32,D21u30,D21u10,D21u12] = D21_coeffs(Rx0,h,D,nu) ;
[D22u20,D22u11,D22u21,D22u31,D22u02,D22u12,D22u22,D22u32,D22u42,D22u13,D22u23,D22u33,D22u24] = D22_coeffs ;
[D2Nu2N,D2Nu2Nm1,D2Nu2Nm2,D2Nu1N,D2Nu3N,D2Nu4N,D2Nu0N,D2Nu3Nm1,D2Nu1Nm1] = D20_coeffs(KxL,RxL,h,D,nu) ;
[D2Nm1u2Nm1,D2Nm1u2Nm2,D2Nm1u2Nm3,D2Nm1u2N,D2Nm1u1Nm1,D2Nm1u3Nm1,D2Nm1u4Nm1,D2Nm1u0Nm1,D2Nm1u3Nm2,D2Nm1u3N,D2Nm1u1N,D2Nm1u1Nm2] = D21_coeffs(RxL,h,D,nu) ;

%%
%%-- Blk31
D0 = D22u02*a0 ;

Blk31               = sparse(diag(D0))   ;
Blk31(1,1)          = D20u00 ;
Blk31(2,2)          = D21u01 ;
Blk31(end,end)      = D2Nu0N ;
Blk31(end-1,end-1)  = D2Nm1u0Nm1 ;

%%
%%-- Blk32
D0 = D22u12*a0 ; D1 = D22u13*a1 ; Dm1 = D22u11*a1 ;

Blk32               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
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

%%
%%-- Blk33
D0 = D22u22*a0 ; D1 = D22u23*a1; D2 = D22u24*a2 ; Dm1 = D22u21*a1; Dm2 = D22u20*a2 ;

Blk33               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
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

%%
%%-- Blk34
D0 = D22u32*a0 ; D1 = D22u33*a1 ; Dm1 = D22u31*a1 ;

Blk34               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
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

%%
%%-- Blk35
D0 = D22u42*a0 ;

Blk35               = sparse(diag(D0))   ;
Blk35(1,1)          = D20u40 ;
Blk35(2,2)          = D21u41 ;
Blk35(end,end)      = D2Nu4N ;
Blk35(end-1,end-1)  = D2Nm1u4Nm1 ;

[D00u00,D00u10,D00u20,D00u01,D00u02,D00u11] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[D01u01,D01u11,D01u21,D01u00,D01u02,D01u03,D01u12,D01u10] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[D02u02,D02u12,D02u22,D02u01,D02u03,D02u04,D02u00,D02u13,D02u11] = D02_coeffs(KLy,RLy,h,D,nu) ;
[D0Nu0N,D0Nu1N,D0Nu2N,D0Nu0Nm1,D0Nu0Nm2,D0Nu1Nm1] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,D0Nm1u1Nm1,D0Nm1u2Nm1,D0Nm1u0N,D0Nm1u0Nm2,D0Nm1u0Nm3,D0Nm1u1Nm2,D0Nm1u1N] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

%%
%%-- BlkMM
D0 = D02u02*a0 ; D1 = D02u03*a1; D2 = D02u04*a2 ; Dm1 = D02u01*a1; Dm2 = D02u00*a2 ;

BlkMM               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
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

%%
%%-- BlkMMm1
D0 = D02u12*a0 ; D1 = D02u13*a1 ; Dm1 = D02u11*a1 ;

BlkMMm1               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
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

%%
%%-- BlkMMm2
D0 = D02u22*a0 ;

BlkMMm2               = sparse(diag(D0))   ;
BlkMMm2(1,1)          = D00u20 ;
BlkMMm2(2,2)          = D01u21 ;
BlkMMm2(end,end)      = D0Nu2N ;
BlkMMm2(end-1,end-1)  = D0Nm1u2Nm1 ;

[D10u10, D10u20, D10u30, D10u00, D10u11, D10u12, D10u21, D10u01] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[D11u11,D11u12,D11u13,D11u10,D11u01,D11u21,D11u31,D11u22,D11u20,D11u00,D11u02] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[D12u12,D12u13,D12u14,D12u11,D12u10,D12u02,D12u22,D12u32,D12u23,D12u21,D12u01,D12u03] = D12_coeffs(RLy,h,D,nu) ;
[D1Nu1N, D1Nu2N, D1Nu3N, D1Nu0N, D1Nu1Nm1, D1Nu1Nm2, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[D1Nm1u1Nm1,D1Nm1u1Nm2,D1Nm1u1Nm3,D1Nm1u1N,D1Nm1u0Nm1,D1Nm1u2Nm1,D1Nm1u3Nm1,D1Nm1u2Nm2,D1Nm1u2N,D1Nm1u0N,D1Nm1u0Nm2] = D11_coeffs(RLy,RxL,h,D,nu) ;

%%
%%-- BlkMm1M
D0 = D12u02*a0 ; D1 = D12u03*a1 ; Dm1 = D12u01*a1 ;

BlkMm1M               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1))  ;
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

%%
%%-- BlkMm1Mm1
D0 = D12u12*a0 ; D1 = D12u13*a1 ; D2 = D12u14*a2 ; Dm1 = D12u11*a1 ; Dm2 = D12u10*a2 ;


BlkMm1Mm1               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1) + diag(D2,2) + diag(Dm2,-2)) ;
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

%%
%%-- BlkMm1Mm2
D0 = D12u22*a0 ; D1 = D12u23*a1 ; Dm1 = D12u21*a1 ;

BlkMm1Mm2               = sparse(diag(D0) + diag(D1,1) + diag(Dm1,-1)) ;
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

%%
%%-- BlkMm1Mm3
D0 = D12u32*a0 ;

BlkMm1Mm3               = sparse(diag(D0))   ;
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

biHarm = biHarm/h^4 ;


%% diagonal case

a0 = ones(Ny+1,1); a1 = ones(Ny,1); a2 = ones(Ny-1,1);

%% dm2Ny  % pad zeros at the end
[~,~,~,~,~,~,D20u00,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D21u01,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,D22u02,~,~,~,~,~,~,~,~] = D22_coeffs ;
[~,~,~,~,~,~,D2Nu0N,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,D2Nm1u0Nm1,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

dm2Ny0 = D22u02*a0;
dm2Ny0([1,2,Ny,Ny+1]) = [D20u00,D21u01,D2Nm1u0Nm1,D2Nu0N];

[~, ~, D10u30, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D11u31,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D12u32,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[~, ~, D1Nu3N, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,D1Nm1u3Nm1,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

dm2Ny1 = D12u32*a0;
dm2Ny1([1,2,Ny,Ny+1]) = [D10u30,D11u31,D1Nm1u3Nm1,D1Nu3N];

[~,~,D00u20,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,D01u21,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,D02u22,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,D0Nu2N,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,D0Nm1u2Nm1,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

dm2Ny2 = D02u22*a0;
dm2Ny2([1,2,Ny,Ny+1]) = [D00u20, D01u21, D0Nm1u2Nm1, D0Nu2N];

Dm2Ny = [repmat(dm2Ny0,Nx-3,1);dm2Ny1;dm2Ny2];

assert(all((diag(Blk31,0) - dm2Ny0) <= eps), "d01 incorrect");
assert(all((diag(BlkMm1Mm3,0) - dm2Ny1) <= eps), "d02 incorrect");
assert(all((diag(BlkMMm2,0) - dm2Ny2) <= eps), "d0Mm1 incorrect");

%% dmNym1 % pad zeros at the end
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,D11u00,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,D12u01,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, D1Nu0Nm1] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,D1Nm1u0Nm2] = D11_coeffs(R0y,RxL,h,D,nu) ;


dmNym10 = D12u01*a1;
dmNym10([1,Ny-1,Ny]) = [D11u00,D1Nm1u0Nm2,D1Nu0Nm1];

[~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,D21u10,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,D22u11,~,~,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
[~,~,~,~,~,~,~,~,D2Nu1Nm1] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~,D2Nm1u1Nm2] = D21_coeffs(RxL,h,D,nu) ;


dmNym11 = D22u11*a1;
dmNym11([1,Ny-1,Ny]) = [D21u10,D2Nm1u1Nm2,D2Nu1Nm1];

[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D01u10] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,D02u11] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,~,~,~,D0Nu1Nm1] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,D0Nm1u1Nm2,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

dmNym1M = D02u11*a1;
dmNym1M([1,Ny-1,Ny]) = [D01u10,D0Nm1u1Nm2,D0Nu1Nm1];

assert(all((diag(Blk21,-1) - dmNym10) <= eps), "d01 incorrect");
assert(all((diag(Blk32,-1) - dmNym11) <= eps), "d02 incorrect");
assert(all((diag(BlkMMm1,-1) - dmNym1M) <= eps), "d0Mm1 incorrect");

DmNym1 = [dmNym10;0;repmat([dmNym11;0],Nx-2,1);dmNym1M];

%% dmNy   % pad zeros at the end

[~, ~, ~, D10u00, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,D11u01,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,D12u02,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, D1Nu0N, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,D1Nm1u0Nm1,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

dmNy0 = D12u02*a0;
dmNy0([1,2,Ny,Ny+1]) = [D10u00,D11u01,D1Nm1u0Nm1,D1Nu0N];

[~,~,~,D20u10,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,D21u11,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,D22u12,~,~,~,~,~,~,~] = D22_coeffs ;
[~,~,~,D2Nu1N,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,D2Nm1u1Nm1,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

dmNy1 = D22u12*a0;
dmNy1([1,2,Ny,Ny+1]) = [D20u10,D21u11,D2Nm1u1Nm1,D2Nu1N];

[~, D10u20, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,D11u21,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D12u22,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[~, D1Nu2N, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,D1Nm1u2Nm1,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

dmNyM1 = D12u22*a0;
dmNyM1([1,2,Ny,Ny+1]) = [D10u20,D11u21,D1Nm1u2Nm1,D1Nu2N];

[~,D00u10,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,D01u11,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,D02u12,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,D0Nu1N,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,D0Nm1u1Nm1,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

dmNyM = D02u12*a0;
dmNyM([1,2,Ny,Ny+1]) = [D00u10,D01u11,D0Nm1u1Nm1,D0Nu1N];

DmNy = [dmNy0;repmat(dmNy1,(Nx-3),1);dmNyM1;dmNyM];

assert(all((diag(Blk21,0) - dmNy0) <= eps), "d01 incorrect");
assert(all((diag(Blk32,0) - dmNy1) <= eps), "d02 incorrect");
assert(all((diag(BlkMMm1,0) - dmNyM) <= eps), "d0Mm1 incorrect");
%assert(all((diag(biHarm,Ny+1) - DmNy) <= eps), "D0 incorrect");


%% dmNyp1 % pad zeros at the
[~, ~, ~, ~, ~, ~, ~, D10u01] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,D11u02] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~,D12u03] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,D1Nm1u0N,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

dmNy10 = D12u03*a1;
dmNy10([1,2,Ny]) = [D10u01,D11u02,D1Nm1u0N];

[~,~,~,~,~,~,~,~,D20u11] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~,D21u12] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,D22u13,~,~,~] = D22_coeffs ;
[~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,D2Nm1u1N,~] = D21_coeffs(RxL,h,D,nu) ;


dmNy11 = D22u13*a1;
dmNy11([1,2,Ny]) = [D20u11,D21u12,D2Nm1u1N];

[~,~,~,~,~,D00u11] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D01u12,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D02u13,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,D0Nm1u1N] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

dmNy1M = D02u13*a1;
dmNy1M([1,2,Ny]) = [D00u11,D01u12,D0Nm1u1N];

assert(all((diag(Blk21,1) - dmNy10) <= eps), "d01 incorrect");
assert(all((diag(Blk32,1) - dmNy11) <= eps), "d02 incorrect");
assert(all((diag(BlkMMm1,1) - dmNy1M) <= eps), "d0Mm1 incorrect");

DmNy1 = [0;dmNy10;repmat([0;dmNy11],Nx-2,1);0;dmNy1M;0];

%% dm2   % pad zeros at the end
[~,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D02u00,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,~,~,~,D0Nu0Nm2,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;


dm20 = D02u00*a2;
dm20([Ny-2,Ny-1])  = [D0Nm1u0Nm3,D0Nu0Nm2];

[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,D12u10,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, ~, D1Nu1Nm2, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

dm21 = D12u10*a2 ;
dm21([Ny-2,Ny-1])  = [D1Nm1u1Nm3,D1Nu1Nm2];

[~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[D22u20,~,~,~,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
[~,~,D2Nu2Nm2,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,D2Nm1u2Nm3,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;


dm22 = D22u20*a2;
dm22([Ny-2,Ny-1])  = [D2Nm1u2Nm3,D2Nu2Nm2];

[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,~,~,~,D12u10,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[~, ~, ~, ~, ~, D1Nu1Nm2, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;


dm2M1 = D12u10*a2 ;
dm2M1([Ny-2,Ny-1]) = [D1Nm1u1Nm3,D1Nu1Nm2];

[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D02u00,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,~,~,D0Nu0Nm2,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

dm2M = D02u00*a2;
dm2M([Ny-2,Ny-1]) = [D0Nm1u0Nm3,D0Nu0Nm2];

assert(all((diag(Blk11,-2) - dm20) <= eps), "d00 incorrect");
assert(all((diag(Blk22,-2) - dm21) <= eps), "d01 incorrect");
assert(all((diag(Blk33,-2) - dm22) <= eps), "d02 incorrect");
assert(all((diag(BlkMm1Mm1,-2) - dm2M1) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMM,-2) - dm2M) <= eps), "d0M incorrect");

Dm2 = [dm20;0;0;dm21;0;0;repmat([dm22;0;0],Nx+1-4,1);dm2M1;0;0;dm2M];

%% dm1   % pad zeros at the end


[~,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,D01u00,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,~,~,D02u01,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,~,~,D0Nu0Nm1,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,D0Nm1u0Nm2,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

dm10 = D02u01*a1;
dm10([1,Ny-1,Ny]) = [D01u00, D0Nm1u0Nm2, D0Nu0Nm1];

[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,D11u10,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,D12u11,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, D1Nu1Nm1, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,D1Nm1u1Nm2,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;


dm11 = D12u11*a1;
dm11([1,Ny-1,Ny]) = [D11u10, D1Nm1u1Nm2, D1Nu1Nm1];

[~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,D21u20,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,D22u21,~,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
[~,D2Nu2Nm1,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,D2Nm1u2Nm2,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;


dm12 = D22u21*a1;
dm12([1,Ny-1,Ny]) = [D21u20, D2Nm1u2Nm2, D2Nu2Nm1];

[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,D11u10,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,~,~,D12u11,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[~, ~, ~, ~, D1Nu1Nm1, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,D1Nm1u1Nm2,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;


dm1M1 = D12u11*a1;
dm1M1([1,Ny-1,Ny]) = [D11u10, D1Nm1u1Nm2, D1Nu1Nm1];

[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,D01u00,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,~,D02u01,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,~,D0Nu0Nm1,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,D0Nm1u0Nm2,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


dm1M = D02u01*a1;
dm1M([1,Ny-1,Ny]) = [D01u00, D0Nm1u0Nm2, D0Nu0Nm1];

assert(all((diag(Blk11,-1) - dm10) <= eps), "d00 incorrect");
assert(all((diag(Blk22,-1) - dm11) <= eps), "d01 incorrect");
assert(all((diag(Blk33,-1) - dm12) <= eps), "d02 incorrect");
assert(all((diag(BlkMm1Mm1,-1) - dm1M1) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMM,-1) - dm1M) <= eps), "d0M incorrect");

Dm1 = [dm10;0;dm11;0;repmat([dm12;0],Nx+1-4,1);dm1M1;0;dm1M];

%% d00
[~,~,~,~,~] = biharmdiag(BCs,h,D,nu);

[D00u00,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[D01u01,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[D02u02,~,~,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[D0Nu0N,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

d00 = D02u02*a0;
d00([1,2,Ny,Ny+1]) = [D00u00, D01u01, D0Nm1u0Nm1, D0Nu0N];

[D10u10, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[D11u11,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[D12u12,~,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[D1Nu1N, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[D1Nm1u1Nm1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;


d01 = D12u12*a0;
d01([1,2,Ny,Ny+1]) = [D10u10, D11u11, D1Nm1u1Nm1, D1Nu1N];

[D20u20,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[D21u21,~,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,D22u22,~,~,~,~,~,~] = D22_coeffs ;
[D2Nu2N,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[D2Nm1u2Nm1,~,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;


d02 = D22u22*a0;
d02([1,2,Ny,Ny+1]) = [D20u20, D21u21, D2Nm1u2Nm1, D2Nu2N];


[D10u10, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[D11u11,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[D12u12,~,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[D1Nu1N, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[D1Nm1u1Nm1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;


d0Mm = D12u12*a0;
d0Mm([1,2,Ny,Ny+1]) = [D10u10, D11u11, D1Nm1u1Nm1, D1Nu1N];

[D00u00,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[D01u01,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[D02u02,~,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[D0Nu0N,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


d0M = D02u02*a0;
d0M([1,2,Ny,Ny+1]) = [D00u00, D01u01, D0Nm1u0Nm1, D0Nu0N];

D0 = [d00;d01;repmat(d02,(Nx+1-4),1);d0Mm;d0M];

assert(all((diag(Blk11,0) - d00) <= eps), "d00 incorrect");
assert(all((diag(Blk22,0) - d01) <= eps), "d01 incorrect");
assert(all((diag(Blk33,0) - d02) <= eps), "d02 incorrect");
assert(all((diag(BlkMm1Mm1,0) - d0Mm) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMM,0) - d0M) <= eps), "d0M incorrect");
% assert(all((diag(biHarm,0) - D0) <= eps), "D0 incorrect");


%% dp1   % pad zeros at the start
[~,~,~,D00u01,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,D01u02,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,~,~,~,D02u03,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,~,~,D0Nm1u0N,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

d10 = D02u03*a1;
d10([1,2,Ny]) = [D00u01, D01u02, D0Nm1u0N];

[~, ~, ~, ~, D10u11, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,D11u12,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,D12u13,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,D1Nm1u1N,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

d11 = D12u13*a1;
d11([1,2,Ny]) = [D10u11, D11u12, D1Nm1u1N];

[~,D20u21,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,D21u22,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,D22u23,~,~] = D22_coeffs ;
[~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,D2Nm1u2N,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

d12 = D22u23*a1;
d12([1,2,Ny]) = [D20u21, D21u22, D2Nm1u2N];

[~, ~, ~, ~, D10u11, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,D11u12,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,D12u13,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,~,~,D1Nm1u1N,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;


d1M1 = D12u13*a1;
d1M1([1,2,Ny]) = [D10u11, D11u12, D1Nm1u1N];

[~,~,~,D00u01,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,D01u02,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,~,~,D02u03,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,~,D0Nm1u0N,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


d1M = D02u03*a1;
d1M([1,2,Ny]) = [D00u01, D01u02, D0Nm1u0N];

D1 = [d10;0;d11;0;repmat([d12;0],Nx+1-4,1);d1M1;0;d1M];

assert(all((diag(Blk11,1) - d10) <= eps), "d00 incorrect");
assert(all((diag(Blk22,1) - d11) <= eps), "d01 incorrect");
assert(all((diag(Blk33,1) - d12) <= eps), "d02 incorrect");
assert(all((diag(BlkMm1Mm1,1) - d1M1) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMM,1) - d1M) <= eps), "d0M incorrect");

%% dp2   % pad zeros at the start
[~,~,~,~,D00u02,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,D01u03,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,D02u04,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

d20 = D02u04*a2;
d20([1,2,Ny-1]) = [D00u02,D01u03,D0Nm1u0Nm3];

[~, ~, ~, ~, ~, D10u12, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,D11u13,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,D12u14,~,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

d21 = D12u14*a2;
d21([1,2,Ny-1]) = [D10u12,D11u13,D1Nm1u1Nm3];

[~,~,D20u22,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,D21u23,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~,~,D22u24] = D22_coeffs ;
[~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,D2Nm1u2Nm3,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

d22 = D22u24*a2;
d22([1,2,Ny-1]) = [D20u22,D21u23,D2Nm1u2Nm3];

[~, ~, ~, ~, ~, D10u12, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,~,D11u13,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,~,D12u14,~,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;


d2M1 = D12u14*a2;
d2M1([1,2,Ny-1]) = [D10u12,D11u13,D1Nm1u1Nm3];

[~,~,~,~,D00u02,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,D01u03,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,D02u04,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

d2M = D02u04*a2;
d2M([1,2,Ny-1]) = [D00u02,D01u03,D0Nm1u0Nm3];

assert(all((diag(Blk11,2) - d20) <= eps), "d00 incorrect");
assert(all((diag(Blk22,2) - d21) <= eps), "d01 incorrect");
assert(all((diag(Blk33,2) - d22) <= eps), "d02 incorrect");
assert(all((diag(BlkMm1Mm1,2) - d2M1) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMM,2) - d2M) <= eps), "d0M incorrect");

D2 = [d20;0;0;d21;0;0;repmat([d22;0;0],Nx+1-4,1);d2M1;0;0;d2M];

%% dpNym1 % pad zeros at the start
[~,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D01u10] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,D02u11] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,~,~,~,~,D0Nu1Nm1] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,D0Nm1u1Nm2,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

dpNym10 = D02u11*a1;
dpNym10([1,Ny-1,Ny]) = [D01u10,D0Nm1u1Nm2,D0Nu1Nm1];

[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,D11u20,D11u00,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,D12u21,D12u01,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, ~, ~, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,D1Nm1u2Nm2,~,~,D1Nm1u0Nm2] = D11_coeffs(R0y,RxL,h,D,nu) ;

dpNym11 = D12u21*a1;
dpNym11([1,Ny-1,Ny]) = [D11u20,D1Nm1u2Nm2,D1Nu2Nm1];

[~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,D21u30,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,D22u31,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
[~,~,~,~,~,~,~,D2Nu3Nm1,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,D2Nm1u3Nm2,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

dpNym12  = D22u31*a1;
dpNym12([1,Ny-1,Ny]) = [D21u30,D2Nm1u3Nm2,D2Nu3Nm1];

[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

dpNym1M = D12u01*a1;
dpNym1M([1,Ny-1,Ny]) = [D11u00,D1Nm1u0Nm2,D1Nu0Nm1];

assert(all((diag(Blk12,-1) - dpNym10) <= eps), "d01 incorrect");
assert(all((diag(Blk23,-1) - dpNym11) <= eps), "d02 incorrect");
assert(all((diag(Blk34,-1) - dpNym12) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMm1M,-1) - dpNym1M) <= eps), "d0Mm1 incorrect");

DNym1 = [0;dpNym10;0;dpNym11;repmat([0;dpNym12],Nx-3,1);0;dpNym1M;0];

%% dpNy   % pad zeros at the start
[~,D00u10,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,D01u11,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,D02u12,~,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,D0Nu1N,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,D0Nm1u1Nm1,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

dpNy0 = D02u12*a0;
dpNy0([1,2,Ny,Ny+1]) = [D00u10,D01u11,D0Nm1u1Nm1,D0Nu1N];

[~, D10u20, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,D11u21,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D12u22,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, D1Nu2N, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,D1Nm1u2Nm1,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;
%%
dpNy1 = D12u22*a0;
dpNy1([1,2,Ny,Ny+1]) = [D10u20,D11u21,D1Nm1u2Nm1,D1Nu2N];


[~,~,~,~,D20u30,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,D21u31,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D22u32,~,~,~,~,~] = D22_coeffs ;
[~,~,~,~,D2Nu3N,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,~,D2Nm1u3Nm1,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

dpNy2 = D22u32*a0;
dpNy2([1,2,Ny,Ny+1]) = [D20u30,D21u31,D2Nm1u3Nm1,D2Nu3N];

[~, ~, ~, D10u00, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,D11u01,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,D12u02,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
[~, ~, ~, D1Nu0N, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,D1Nm1u0Nm1,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

dpNyM = D12u02*a0;
dpNyM([1,2,Ny,Ny+1]) = [D10u00,D11u01,D1Nm1u0Nm1,D1Nu0N];

DNy = [dpNy0;dpNy1;repmat(dpNy2,(Nx-3),1);dpNyM];

assert(all((diag(Blk12,0) - dpNy0) <= eps), "d01 incorrect");
assert(all((diag(Blk23,0) - dpNy1) <= eps), "d02 incorrect");
assert(all((diag(Blk34,0) - dpNy2) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMm1M,0) - dpNyM) <= eps), "d0Mm1 incorrect");

%% dpNyp1 % pad zeros at the start
[~,~,~,~,~,D00u11] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D01u12,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D02u13,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,D0Nm1u1N] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

dpNy10 = D02u13*a1 ;
dpNy10([1,2,Ny]) = [D00u11, D01u12, D0Nm1u1N];

[~, ~, ~, ~, ~, ~, D10u21, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D11u22,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,D12u23,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,D1Nm1u2N,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;


dpNy11 = D12u23*a1 ;
dpNy11([1,2,Ny]) = [D10u21, D11u22, D1Nm1u2N];

[~,~,~,~,~,~,~,D20u31,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,D21u32,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~,D22u33,~] = D22_coeffs ;
[~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,D2Nm1u3N,~,~] = D21_coeffs(RxL,h,D,nu) ;

dpNy12 = D22u33*a1;
dpNy12([1,2,Ny]) = [D20u31, D21u32, D2Nm1u3N];

[~, ~, ~, ~, ~, ~, ~, D10u01] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,D11u02] = D11_coeffs(RLy,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,~,~,D12u03] = D12_coeffs(RLy,h,D,nu) ;
[~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,~,~,~,D1Nm1u0N,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

dpNy1M = D12u03*a1;
dpNy1M([1,2,Ny]) = [D10u01, D11u02, D1Nm1u0N];

assert(all((diag(Blk12,1) - dpNy10) <= eps), "d01 incorrect");
assert(all((diag(Blk23,1) - dpNy11) <= eps), "d02 incorrect");
assert(all((diag(Blk34,1) - dpNy12) <= eps), "d0Mm1 incorrect");
assert(all((diag(BlkMm1M,1) - dpNy1M) <= eps), "d0Mm1 incorrect");

DNy1 = [dpNy10;0;dpNy11;0;repmat([dpNy12;0],Nx-3,1);dpNy1M];

%% dp2Ny  % pad zeros at the start
[~,~,D00u20,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[~,~,D01u21,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[~,~,D02u22,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
[~,~,D0Nu2N,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[~,~,D0Nm1u2Nm1,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

dp2Ny0 = D02u22*a0;
dp2Ny0([1,2,Ny,Ny+1]) = [D00u20,D01u21,D0Nm1u2Nm1,D0Nu2N];

[~, ~, D10u30, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D11u31,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,D12u32,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
[~, ~, D1Nu3N, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,D1Nm1u3Nm1,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

dp2Ny1 = D12u32*a0;
dp2Ny1([1,2,Ny,Ny+1]) = [D10u30,D11u31,D1Nm1u3Nm1,D1Nu3N];

[~,~,~,~,~,D20u40,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[~,~,~,~,~,~,D21u41,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,~,~,D22u42,~,~,~,~] = D22_coeffs ;
[~,~,~,~,~,D2Nu4N,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
[~,~,~,~,~,~,D2Nm1u4Nm1,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

dp2Ny2 = D22u42*a0;
dp2Ny2([1,2,Ny,Ny+1]) = [D20u40 ,D21u41 ,D2Nm1u4Nm1 ,D2Nu4N];

D2Ny = [dp2Ny0;dp2Ny1;repmat(dp2Ny2,Nx-3,1)];

assert(all((diag(Blk13,0) - dp2Ny0) <= eps), "d01 incorrect");
assert(all((diag(Blk24,0) - dp2Ny1) <= eps), "d02 incorrect");
assert(all((diag(Blk35,0) - dp2Ny2) <= eps), "d0Mm1 incorrect");

%% Zero Padding
Dm2Ny = [Dm2Ny; zeros((Nx+1)*2,1)];

DmNym1 = [DmNym1; zeros((Nx+2),1)];
DmNy = [DmNy; zeros((Nx + 1),1)];
DmNy1 = [DmNy1; zeros((Nx),1) ];

Dm2 = [Dm2;0;0];
Dm1 = [Dm1;0];
D1 = [0;D1];
D2 = [0;0;D2];

DNym1 = [zeros((Nx),1);DNym1];
DNy = [zeros((Nx+1),1);DNy]; 
DNy1 = [zeros((Nx+2),1);DNy1];

D2Ny = [zeros((Nx+1)*2,1);D2Ny];
%% diag biharmonic

BHdiags = [Dm2Ny,...
    DmNym1, DmNy, DmNy1,...
    Dm2, Dm1, D0, D1, D2,...
    DNym1, DNy, DNy1,...
    D2Ny];

dn = [-(2*(Nx+1)),-(Nx+2),-(Nx+1),-(Nx), (-2:2), (Nx),(Nx+1),(Nx+2), 2*(Nx+1)];

BH = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
BH = spdiags(BHdiags, dn, BH)/h^4;

spy(biHarm - BH);

end

