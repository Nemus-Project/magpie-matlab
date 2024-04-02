function [d0,d1,dc,dM1,dM] = biharmdiag(BCs,h,D,nu,dia)
%BIHARMDIAG Summary of this function goes here
%   Detailed explanation goes here

pack_BCs = num2cell(BCs);
[K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = pack_BCs{:};

switch dia
    case '-2ny'
        [~,~,~,~,~,~,D20u00,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D21u01,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,D22u02,~,~,~,~,~,~,~,~] = D22_coeffs ;
        [~,~,~,~,~,~,~,D2Nm1u0Nm1,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;
        [~,~,~,~,~,~,D2Nu0N,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;

        [~, ~, D10u30, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D11u31,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D12u32,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~,~,~,~,~,~,D1Nm1u3Nm1,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;
        [~, ~, D1Nu3N, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;        

        [~,~,D00u20,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,D01u21,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,D02u22,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,D0Nm1u2Nm1,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;
        [~,~,D0Nu2N,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        
        d0 = [NaN; NaN; NaN; NaN; NaN];
        d1 = [NaN; NaN; NaN; NaN; NaN];
        dc = [D20u00; D21u01; D22u02; D2Nm1u0Nm1; D2Nu0N];
        dM1 = [D10u30; D11u31; D12u32; D1Nm1u3Nm1; D1Nu3N];
        dM = [D00u20; D01u21; D02u22; D0Nm1u2Nm1; D0Nu2N];
        
    case '-ny-1'

        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,D11u00,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,D12u01,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, D1Nu0Nm1] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,D1Nm1u0Nm2] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,D21u10,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,D22u11,~,~,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
        [~,~,~,~,~,~,~,~,D2Nu1Nm1] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~,D2Nm1u1Nm2] = D21_coeffs(RxL,h,D,nu) ;

        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D01u10] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D02u11] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,~,~,~,D0Nu1Nm1] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,D0Nm1u1Nm2,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

        
        d0 = [D11u00; D12u01; D1Nm1u0Nm2; D1Nu0Nm1;];
        d1 = [NaN; NaN; NaN; NaN;];
        dc = [D21u10; D22u11; D2Nm1u1Nm2; D2Nu1Nm1];
        dM1 = [NaN; NaN; NaN; NaN;];
        dM = [D01u10; D02u11; D0Nm1u1Nm2; D0Nu1Nm1];
        
        

    case '-ny'

        [~, ~, ~, D10u00, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,D11u01,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,D12u02,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, D1Nu0N, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,D1Nm1u0Nm1,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,D20u10,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,D21u11,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,D22u12,~,~,~,~,~,~,~] = D22_coeffs ;
        [~,~,~,D2Nu1N,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,~,D2Nm1u1Nm1,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~, D10u20, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,D11u21,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D12u22,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~, D1Nu2N, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,D1Nm1u2Nm1,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

        [~,D00u10,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,D01u11,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,D02u12,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,D0Nu1N,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,D0Nm1u1Nm1,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


    case '-ny+1'

        [~, ~, ~, ~, ~, ~, ~, D10u01] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,D11u02] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~,D12u03] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,D1Nm1u0N,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,~,~,~,~,D20u11] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~,D21u12] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,D22u13,~,~,~] = D22_coeffs ;
        [~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,D2Nm1u1N,~] = D21_coeffs(RxL,h,D,nu) ;

        [~,~,~,~,~,D00u11] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D01u12,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D02u13,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,D0Nm1u1N] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


    case '-2'

        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D02u00,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,~,~,~,D0Nu0Nm2,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,D12u10,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, ~, D1Nu1Nm2, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [D22u20,~,~,~,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
        [~,~,D2Nu2Nm2,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,D2Nm1u2Nm3,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,D12u10,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~, ~, ~, ~, ~, D1Nu1Nm2, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
        [~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D02u00,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,~,~,D0Nu0Nm2,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


    case '-1'

        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,D01u00,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,~,~,D02u01,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,~,~,D0Nu0Nm1,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,D0Nm1u0Nm2,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,D11u10,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,D12u11,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, D1Nu1Nm1, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,D1Nm1u1Nm2,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,D21u20,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,D22u21,~,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
        [~,D2Nu2Nm1,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,D2Nm1u2Nm2,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,D11u10,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,D12u11,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~, ~, ~, ~, D1Nu1Nm1, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
        [~,D1Nm1u1Nm2,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,D01u00,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,~,D02u01,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,~,D0Nu0Nm1,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,D0Nm1u0Nm2,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

    case '0'

        [D00u00,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [D01u01,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [D02u02,~,~,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu);
        [D0Nu0N,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [D0Nm1u0Nm1,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [D10u10,~,~,~,~,~,~,~]           = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [D11u11,~,~,~,~,~,~,~,~,~,~]     = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [D12u12,~,~,~,~,~,~,~,~,~,~,~]   = D12_coeffs(R0y,h,D,nu) ;
        [D1Nu1N,~,~,~,~,~,~,~]           = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [D1Nm1u1Nm1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [D20u20,~,~,~,~,~,~,~,~]           = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [D21u21,~,~,~,~,~,~,~,~,~,~,~]     = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D22u22,~,~,~,~,~,~]   = D22_coeffs ;
        [D2Nu2N,~,~,~,~,~,~,~,~]           = D20_coeffs(KxL,RxL,h,D,nu) ;
        [D2Nm1u2Nm1,~,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [D10u10,~,~,~,~,~,~,~]    = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [D11u11,~,~,~,~,~,~,~,~,~,~]     = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [D12u12,~,~,~,~,~,~,~,~,~,~,~]   = D12_coeffs(RLy,h,D,nu) ;
        [D1Nm1u1Nm1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;
        [D1Nu1N,~,~,~,~,~,~,~]    = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;

        [D00u00,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [D01u01,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [D02u02,~,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [D0Nu0N,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [D0Nm1u0Nm1,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

        d0 = [D00u00; D01u01; D02u02; D0Nm1u0Nm1; D0Nu0N];
        d1 = [D10u10; D11u11; D12u12; D1Nm1u1Nm1; D1Nu1N];
        dc = [D20u20; D21u21; D22u22; D2Nm1u2Nm1; D2Nu2N];
        dM1 = [D10u10; D11u11; D12u12; D1Nm1u1Nm1; D1Nu1N];
        dM = [D00u00; D01u01; D02u02; D0Nm1u0Nm1; D0Nu0N];



    case '1'
        [~,~,~,D00u01,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,D01u02,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,~,~,~,D02u03,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,D0Nm1u0N,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, ~, ~, ~, D10u11, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,D11u12,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,D12u13,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,D1Nm1u1N,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,D20u21,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,D21u22,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,D22u23,~,~] = D22_coeffs ;
        [~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,D2Nm1u2N,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~, ~, ~, ~, D10u11, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,D11u12,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,D12u13,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,D1Nm1u1N,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

        [~,~,~,D00u01,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,D01u02,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,~,~,D02u03,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,D0Nm1u0N,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


    case '2'

        [~,~,~,~,D00u02,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,D01u03,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,D02u04,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, D10u12, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,D11u13,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,D12u14,~,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,D20u22,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,D21u23,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~,~,D22u24] = D22_coeffs ;
        [~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,D2Nm1u2Nm3,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, D10u12, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,D11u13,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,D12u14,~,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
        [~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

        [~,~,~,~,D00u02,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,D01u03,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,D02u04,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;


    case 'ny-1'

        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D01u10] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D02u11] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,~,~,~,~,D0Nu1Nm1] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,D0Nm1u1Nm2,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D11u20,D11u00,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,D12u21,D12u01,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, D1Nu2Nm1, D1Nu0Nm1] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,D1Nm1u2Nm2,~,~,D1Nm1u0Nm2] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,D21u30,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,D22u31,~,~,~,~,~,~,~,~,~] = D22_coeffs ;
        [~,~,~,~,~,~,~,D2Nu3Nm1,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D2Nm1u3Nm2,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;



    case 'ny'

        [~,D00u10,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,D01u11,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,D02u12,~,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,D0Nu1N,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,D0Nm1u1Nm1,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, D10u20, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,D11u21,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D12u22,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, D1Nu2N, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,D1Nm1u2Nm1,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,D20u30,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,D21u31,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D22u32,~,~,~,~,~] = D22_coeffs ;
        [~,~,~,~,D2Nu3N,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,D2Nm1u3Nm1,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~, ~, ~, D10u00, ~, ~, ~, ~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,D11u01,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,D12u02,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~, ~, ~, D1Nu0N, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,D1Nm1u0Nm1,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;


    case 'ny+1'

        [~,~,~,~,~,D00u11] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D01u12,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D02u13,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,D0Nm1u1N] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, ~, D10u21, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D11u22,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D12u23,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D1Nm1u2N,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,~,~,~,D20u31,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D21u32,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~,D22u33,~] = D22_coeffs ;
        [~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,D2Nm1u3N,~,~] = D21_coeffs(RxL,h,D,nu) ;

        [~, ~, ~, ~, ~, ~, ~, D10u01] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,D11u02] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,~,D12u03] = D12_coeffs(RLy,h,D,nu) ;
        [~, ~, ~, ~, ~, ~, ~, ~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,D1Nm1u0N,~] = D11_coeffs(RLy,RxL,h,D,nu) ;


    case '2ny'

        [~,~,D00u20,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,D01u21,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
        [~,~,D02u22,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu) ;
        [~,~,D0Nu2N,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
        [~,~,D0Nm1u2Nm1,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

        [~, ~, D10u30, ~, ~, ~, ~, ~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D11u31,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,D12u32,~,~,~,~] = D12_coeffs(R0y,h,D,nu) ;
        [~, ~, D1Nu3N, ~, ~, ~, ~, ~] = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,D1Nm1u3Nm1,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

        [~,~,~,~,~,D20u40,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,D21u41,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,D22u42,~,~,~,~] = D22_coeffs ;
        [~,~,~,~,~,D2Nu4N,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu) ;
        [~,~,~,~,~,~,D2Nm1u4Nm1,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;



end

D0 = [d0,d1,dc,dM1,dM];

end

