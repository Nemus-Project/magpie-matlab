function [d0,d1,dc,dM1,dM] = biharmdiag(BCs,h,D,nu,dia)
%BIHARMDIAG Summary of this function goes here
%   Detailed explanation goes here

pack_BCs = num2cell(BCs);
[K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = pack_BCs{:};

d0 = NaN;
d1 = NaN;
dc = NaN;
dM1 = NaN;
dM = NaN;


switch dia
    case '-2ny'
        % if(method)
        %     D20u00 = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + Rx0*h)) - (4*D*nu)/(2*D + Rx0*h) + 1);
        %     D21u01 = 2;
        %     D22u02 = 1
        %     D2Nm1u0Nm1 = 2
        %     D2Nu0N = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + RxL*h)) - (4*D*nu)/(2*D + RxL*h) + 1);

        %     D10u30 = ((16*D^4*nu - 8*D^4*nu^2 - 16*D^4*nu^3 + 8*D^4*nu^4 +
        %     8*D^3*RLy*h*nu + 8*D^3*Rx0*h*nu - 4*D^3*RLy*h*nu^2 -
        %     4*D^3*Rx0*h*nu^2 - 2*D^2*RLy*Rx0*h^2*nu^2 +
        %     4*D^2*RLy*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 +
        %     2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2)) - (4*D*nu)/(2*D + Rx0*h) +
        %     1);

        %     D11u31 = 1
        %     D12u32 = 1
        %     D1Nm1u3Nm1 = 1
        %     D1Nu3N = ((16*D^4*nu - 8*D^4*nu^2 - 16*D^4*nu^3 + 8*D^4*nu^4 + 8*D^3*RLy*h*nu + 8*D^3*RxL*h*nu - 4*D^3*RLy*h*nu^2 - 4*D^3*RxL*h*nu^2 - 2*D^2*RLy*RxL*h^2*nu^2 + 4*D^2*RLy*RxL*h^2*nu)/(D*(2*D + RxL*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2)) - (4*D*nu)/(2*D + RxL*h) + 1);

        % else
        [~,~,~,~,~,~,D20u00,~,~]          = D20_coeffs(Kx0,Rx0,h,D,nu,6);
        [~,~,~,~,~,~,~,D21u01,~,~,~,~]    = D21_coeffs(Rx0,h,D,nu,8);
        [~,~,~,~,D22u02,~,~,~,~,~,~,~,~]  = D22_coeffs;
        [~,~,~,~,~,~,~,D2Nm1u0Nm1,~,~,~,~]= D21_coeffs(RxL,h,D,nu,8);
        [~,~,~,~,~,~,D2Nu0N,~,~]          = D20_coeffs(KxL,RxL,h,D,nu,6);

        [~,~, D10u30,~,~,~,~,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu,2);
        [~,~,~,~,~,~,D11u31,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu,6);
        [~,~,~,~,~,~,~,D12u32,~,~,~,~] = D12_coeffs(RLy,h,D,nu,7);
        [~,~,~,~,~,~,D1Nm1u3Nm1,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu,6);
        [~,~, D1Nu3N,~,~,~,~,~] = D10_coeffs(RLy,KxL,RxL,h,D,nu,2);

        [~,~,D00u20,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,2);
        [~,~,D01u21,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu,2);
        [~,~,D02u22,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu,2);
        [~,~,D0Nm1u2Nm1,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu,2);
        [~,~,D0Nu2N,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,2);

        d0 = [NaN; NaN; NaN; NaN; NaN];
        d1 = [NaN; NaN; NaN; NaN; NaN];
        dc = [D20u00; D21u01; D22u02; D2Nm1u0Nm1; D2Nu0N];
        dM1 = [D10u30; D11u31; D12u32; D1Nm1u3Nm1; D1Nu3N];
        dM = [D00u20; D01u21; D02u22; D0Nm1u2Nm1; D0Nu2N];

    case '-ny-1'

        [~,~,~,~,~,~,~,~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,~,D11u00,~] = D11_coeffs(R0y,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,~,~,D12u01,~] = D12_coeffs(R0y,h,D,nu);
        [~,~,~,~,~,~,~,~,~,~,D1Nm1u0Nm2] = D11_coeffs(R0y,RxL,h,D,nu);
        [~,~,~,~,~,~,~, D1Nu0Nm1] = D10_coeffs(R0y,KxL,RxL,h,D,nu,7);

        [~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,~,~,D21u10,~] = D21_coeffs(Rx0,h,D,nu);
        [~,D22u11,~,~,~,~,~,~,~,~,~,~,~] = D22_coeffs;
        [~,~,~,~,~,~,~,~,~,~,~,D2Nm1u1Nm2] = D21_coeffs(RxL,h,D,nu);
        [~,~,~,~,~,~,~,~,D2Nu1Nm1] = D20_coeffs(KxL,RxL,h,D,nu,8);

        [~,~,~,~,~,~,~,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,DN11u20,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,DN12u21,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~,~,~,~,~,~,~,DN1Nm1u2Nm2,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;
        [~,~,~,~,~,~, DN1Nu2Nm1,~] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;

        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,D01u10] = D01_coeffs(KLy,RLy,Rx0,h,D,nu,7);
        [~,~,~,~,~,~,~,~,D02u11] = D02_coeffs(KLy,RLy,h,D,nu,8);
        [~,~,~,~,~,~,D0Nm1u1Nm2,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu,6);
        [~,~,~,~,~,D0Nu1Nm1] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,5);

        d0 = [NaN; NaN; NaN; NaN;];
        d1 = [D11u00; D12u01; D1Nm1u0Nm2; D1Nu0Nm1;];
        dc = [D21u10; D22u11; D2Nm1u1Nm2; D2Nu1Nm1];
        dM1 = [DN11u20; DN12u21; DN1Nm1u2Nm2; DN1Nu2Nm1];
        dM = [D01u10; D02u11; D0Nm1u1Nm2; D0Nu1Nm1];



    case '-ny'

        [~,~,~, D10u00,~,~,~,~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,3);
        [~,~,~,~,D11u01,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu,4);
        [~,~,~,~,~,D12u02,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu,5);
        [~,~,~,~,D1Nm1u0Nm1,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu,4);
        [~,~,~, D1Nu0N,~,~,~,~] = D10_coeffs(R0y,KxL,RxL,h,D,nu,3);

        [~,~,~,D20u10,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu,3);
        [~,~,~,~,D21u11,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu,4);
        [~,~,~,~,~,D22u12,~,~,~,~,~,~,~] = D22_coeffs;
        [~,~,~,~,D2Nm1u1Nm1,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu,4);
        [~,~,~,D2Nu1N,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu,3);

        [~, D10u20,~,~,~,~,~,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu,1);
        [~,~,~,~,~,D11u21,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu,6);
        [~,~,~,~,~,~,D12u22,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu,7);
        [~,~,~,~,~,D1Nm1u2Nm1,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu,6);
        [~, D1Nu2N,~,~,~,~,~,~] = D10_coeffs(RLy,KxL,RxL,h,D,nu,1);

        [~,D00u10,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,1);
        [~,D01u11,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu,1);
        [~,D02u12,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu,1);
        [~,D0Nm1u1Nm1,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu,1);
        [~,D0Nu1N,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,1);

        d0 = [NaN; NaN; NaN; NaN;];
        d1 = [D10u00;D11u01;D12u02;D1Nm1u0Nm1;D1Nu0N];
        dc = [D20u10;D21u11;D22u12;D2Nm1u1Nm1;D2Nu1N];
        dM1 = [D10u20;D11u21;D12u22;D1Nm1u2Nm1;D1Nu2N];
        dM = [D00u10;D01u11;D02u12;D0Nm1u1Nm1;D0Nu1N];

    case '-ny+1'

        [~,~,~,~,~,~,~, D10u01] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,7);
        [~,~,~,~,~,~,~,~,~,~,D11u02] = D11_coeffs(R0y,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,~,~,~,D12u03] = D12_coeffs(R0y,h,D,nu);
        [~,~,~,~,~,~,~,~,~,D1Nm1u0N,~] = D11_coeffs(R0y,RxL,h,D,nu);

        [~,~,~,~,~,~,~,~,D20u11] = D20_coeffs(Kx0,Rx0,h,D,nu,8);
        [~,~,~,~,~,~,~,~,~,~,~,D21u12] = D21_coeffs(Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,~,D22u13,~,~,~] = D22_coeffs;
        [~,~,~,~,~,~,~,~,~,~,D2Nm1u1N,~] = D21_coeffs(RxL,h,D,nu);

        [~,~,~,~,~,~,DN10u21,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,DN11u22,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,DN12u23,~,~,~] = D12_coeffs(RLy,h,D,nu) ;
        [~,~,~,~,~,~,~,~,DN1Nm1u2N,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;

        [~,~,~,~,~,D00u11] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,5);
        [~,~,~,~,~,~,D01u12,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu,6);
        [~,~,~,~,~,~,~,D02u13,~] = D02_coeffs(KLy,RLy,h,D,nu,7);
        [~,~,~,~,~,~,~,D0Nm1u1N] = D01_coeffs(KLy,RLy,RxL,h,D,nu,7);

        d0 = [NaN; NaN; NaN; NaN;];
        
        d1 = [D10u01;D11u02;D12u03;D1Nm1u0N];        
        dc = [D20u11;D21u12;D22u13;D2Nm1u1N];
        dM1 = [DN10u21;DN11u22;DN12u23;DN1Nm1u2N];
        dM = [D00u11;D01u12;D02u13;D0Nm1u1N];

    case '-2'

        [~,~,~,~,~,~,D02u00,~,~] = D02_coeffs(K0y,R0y,h,D,nu,6);
        [~,~,~,~,~,D0Nm1u0Nm3,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu,5);
        [~,~,~,~,D0Nu0Nm2,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,4);

        [~,~,~,~,D12u10,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu,4);
        [~,~,D1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu,2);
        [~,~,~,~,~, D1Nu1Nm2,~,~] = D10_coeffs(R0y,KxL,RxL,h,D,nu,5);

        [D22u20,~,~,~,~,~,~,~,~,~,~,~,~] = D22_coeffs;
        [~,~,D2Nm1u2Nm3,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu,2);
        [~,~,D2Nu2Nm2,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu,2);

        [~,~,~,~,DN12u10,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu,4);
        [~,~,DN1Nm1u1Nm3,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu,2);
        [~,~,~,~,~,DN1Nu1Nm2,~,~] = D10_coeffs(RLy,KxL,RxL,h,D,nu,5);

        [~,~,~,~,~,~,DN2u00,~,~] = D02_coeffs(KLy,RLy,h,D,nu,6);
        [~,~,~,~,~,DNNm1u0Nm3,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu,5);
        [~,~,~,~,DNNu0Nm2,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,4);

        d0 = [D02u00; D0Nm1u0Nm3; D0Nu0Nm2];
        d1 = [D12u10;D1Nm1u1Nm3;D1Nu1Nm2];
        dc = [D22u20;D2Nm1u2Nm3;D2Nu2Nm2];
        dM1 = [DN12u10;DN1Nm1u1Nm3;DN1Nu1Nm2];
        dM = [DN2u00;DNNm1u0Nm3;DNNu0Nm2];

    case '-1'

        [~,~,~,~,~,~]                   = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu);
        [~,~,~,D01u00,~,~,~,~]          = D01_coeffs(K0y,R0y,Rx0,h,D,nu,3);
        [~,~,~,D02u01,~,~,~,~,~]        = D02_coeffs(K0y,R0y,h,D,nu,3);
        [~,~,~,~,D0Nm1u0Nm2,~,~,~]      = D01_coeffs(K0y,R0y,RxL,h,D,nu,4);
        [~,~,~,D0Nu0Nm1,~,~]            = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,3);

        [~,~,~,~,~,~,~,~]           = D10_coeffs(R0y,Kx0,Rx0,h,D,nu);
        [~,~,~,D11u10,~,~,~,~,~,~,~]       = D11_coeffs(R0y,Rx0,h,D,nu,3);
        [~,~,~,D12u11,~,~,~,~,~,~,~,~]     = D12_coeffs(R0y,h,D,nu,3);
        [~,D1Nm1u1Nm2,~,~,~,~,~,~,~,~,~]   = D11_coeffs(R0y,RxL,h,D,nu,1);
        [~,~,~,~, D1Nu1Nm1,~,~,~]    = D10_coeffs(R0y,KxL,RxL,h,D,nu,4);

        [~,~,~,~,~,~,~,~,~]                = D20_coeffs(Kx0,Rx0,h,D,nu);
        [~,~,~,D21u20,~,~,~,~,~,~,~,~]     = D21_coeffs(Rx0,h,D,nu,3);
        [~,~,D22u21,~,~,~,~,~,~,~,~,~,~]   = D22_coeffs;
        [~,D2Nm1u2Nm2,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu,1);
        [~,D2Nu2Nm1,~,~,~,~,~,~,~]         = D20_coeffs(KxL,RxL,h,D,nu,1);

        [~,~,~,~,~,~,~,~]         = D10_coeffs(RLy,Kx0,Rx0,h,D,nu);
        [~,~,~,DN11u10,~,~,~,~,~,~,~]     = D11_coeffs(RLy,Rx0,h,D,nu,3);
        [~,~,~,DN12u11,~,~,~,~,~,~,~,~]   = D12_coeffs(RLy,h,D,nu,3);
        [~,DN1Nm1u1Nm2,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu,1);
        [~,~,~,~, DN1Nu1Nm1,~,~,~]  = D10_coeffs(RLy,KxL,RxL,h,D,nu,4);

        [~,~,~,~,~,~]               = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu);
        [~,~,~,DN1u00,~,~,~,~]      = D01_coeffs(KLy,RLy,Rx0,h,D,nu,3);
        [~,~,~,DN2u01,~,~,~,~,~]    = D02_coeffs(KLy,RLy,h,D,nu,3);
        [~,~,~,~,DNNm1u0Nm2,~,~,~]  = D01_coeffs(KLy,RLy,RxL,h,D,nu,4);
        [~,~,~,DNNu0Nm1,~,~]        = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,3);

        d0  = [D01u00;D02u01;D0Nm1u0Nm2;D0Nu0Nm1];
        d1  = [D11u10;D12u11;D1Nm1u1Nm2;D1Nu1Nm1];
        dc  = [D21u20;D22u21;D2Nm1u2Nm2;D2Nu2Nm1];
        dM1 = [DN11u10;DN12u11;DN1Nm1u1Nm2;DN1Nu1Nm1];
        dM  = [DN1u00;DN2u01;DNNm1u0Nm2;DNNu0Nm1];

    case '0'

        % if(method)
        % B0u0 = ((2*K0y*Rx0*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*Rx0*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 12*Rx0*D^2*R0y*h^2*nu^2 + 24*Rx0*D^2*R0y*h^2*nu + 24*Rx0*D^2*R0y*h^2 + 16*D^3*R0y*h*nu^3 - 56*D^3*R0y*h*nu^2 + 48*D^3*R0y*h + 8*K0y*Rx0*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 24*Rx0*D^3*h*nu^2 + 48*Rx0*D^3*h*nu + 48*Rx0*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*(- 8*D^2*nu^2 + 4*Rx0*h*D*nu + 8*D^2 + 4*Rx0*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (4*D*nu)/(2*D + R0y*h) - (4*D*nu)/(2*D + Rx0*h) - (2*(8*D^2*nu + 2*D*R0y*h*nu + 2*D*Rx0*h*nu))/((2*D + R0y*h)*(2*D + Rx0*h)) - (8*(- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (2*Kx0*R0y*Rx0^2*h^6 + 4*Kx0*D*Rx0^2*h^5 + 8*Kx0*R0y*D*Rx0*h^5 - 8*Kx0*D^2*Rx0*h^4*nu^2 + 16*Kx0*D^2*Rx0*h^4 - 12*R0y*D^2*Rx0*h^2*nu^2 + 24*R0y*D^2*Rx0*h^2*nu + 24*R0y*D^2*Rx0*h^2 + 16*D^3*Rx0*h*nu^3 - 56*D^3*Rx0*h*nu^2 + 48*D^3*Rx0*h + 8*Kx0*R0y*D^2*h^4 - 16*Kx0*D^3*h^3*nu^2 + 16*Kx0*D^3*h^3 - 24*R0y*D^3*h*nu^2 + 48*R0y*D^3*h*nu + 48*R0y*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + 20);;
        % B0u1  = ((4*D^2*nu^2 - 4*D^2 - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) - (4*D*nu)/(2*D + R0y*h) + (2*K0y*Rx0*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*Rx0*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 14*Rx0*D^2*R0y*h^2*nu^2 + 28*Rx0*D^2*R0y*h^2*nu + 24*Rx0*D^2*R0y*h^2 - 20*D^3*R0y*h*nu^2 + 40*D^3*R0y*h*nu + 48*D^3*R0y*h + 8*K0y*Rx0*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 28*Rx0*D^3*h*nu^2 + 56*Rx0*D^3*h*nu + 48*Rx0*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*Rx0*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20);;
        % B0u2 = ((2*K0y*R0y*h^4 + 4*K0y*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + R0y*h)) - (8*D*nu)/(2*D + R0y*h) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) + 20);;
        % B0uN1  = ((4*D^2*nu^2 - 4*D^2 - 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) - (4*D*nu)/(2*D + R0y*h) + (2*K0y*RxL*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*RxL*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 14*RxL*D^2*R0y*h^2*nu^2 + 28*RxL*D^2*R0y*h^2*nu + 24*RxL*D^2*R0y*h^2 - 20*D^3*R0y*h*nu^2 + 40*D^3*R0y*h*nu + 48*D^3*R0y*h + 8*K0y*RxL*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 28*RxL*D^3*h*nu^2 + 56*RxL*D^3*h*nu + 48*RxL*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)) - (8*D*RxL*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) + 20);;
        % B0uN = ((2*K0y*RxL*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*RxL*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 12*RxL*D^2*R0y*h^2*nu^2 + 24*RxL*D^2*R0y*h^2*nu + 24*RxL*D^2*R0y*h^2 + 16*D^3*R0y*h*nu^3 - 56*D^3*R0y*h*nu^2 + 48*D^3*R0y*h + 8*K0y*RxL*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 24*RxL*D^3*h*nu^2 + 48*RxL*D^3*h*nu + 48*RxL*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)) - (8*(- 8*D^2*nu^2 + 4*RxL*h*D*nu + 8*D^2 + 4*RxL*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) - (4*D*nu)/(2*D + R0y*h) - (4*D*nu)/(2*D + RxL*h) - (2*(8*D^2*nu + 2*D*R0y*h*nu + 2*D*RxL*h*nu))/((2*D + R0y*h)*(2*D + RxL*h)) - (8*(- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) + (2*KxL*R0y*RxL^2*h^6 + 4*KxL*D*RxL^2*h^5 + 8*KxL*R0y*D*RxL*h^5 - 8*KxL*D^2*RxL*h^4*nu^2 + 16*KxL*D^2*RxL*h^4 - 12*R0y*D^2*RxL*h^2*nu^2 + 24*R0y*D^2*RxL*h^2*nu + 24*R0y*D^2*RxL*h^2 + 16*D^3*RxL*h*nu^3 - 56*D^3*RxL*h*nu^2 + 48*D^3*RxL*h + 8*KxL*R0y*D^2*h^4 - 16*KxL*D^3*h^3*nu^2 + 16*KxL*D^3*h^3 - 24*R0y*D^3*h*nu^2 + 48*R0y*D^3*h*nu + 48*R0y*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + RxL*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)) + 20);;

        % B1u0  = ((4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) - (4*D*nu)/(2*D + Rx0*h) + (2*Kx0*R0y*Rx0^2*h^6 + 4*Kx0*D*Rx0^2*h^5 + 8*Kx0*R0y*D*Rx0*h^5 - 8*Kx0*D^2*Rx0*h^4*nu^2 + 16*Kx0*D^2*Rx0*h^4 - 14*R0y*D^2*Rx0*h^2*nu^2 + 28*R0y*D^2*Rx0*h^2*nu + 24*R0y*D^2*Rx0*h^2 - 20*D^3*Rx0*h*nu^2 + 40*D^3*Rx0*h*nu + 48*D^3*Rx0*h + 8*Kx0*R0y*D^2*h^4 - 16*Kx0*D^3*h^3*nu^2 + 16*Kx0*D^3*h^3 - 28*R0y*D^3*h*nu^2 + 56*R0y*D^3*h*nu + 48*R0y*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20);;
        % B1u1  = (20 - (2*D - Rx0*h)/(2*D + Rx0*h) - (2*D - R0y*h)/(2*D + R0y*h));
        % B1u2  = (20 - (2*D - R0y*h)/(2*D + R0y*h));
        % B1uN1 = ((4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) - (4*D*nu)/(2*D + Rx0*h) + (2*RxL*R0y*Rx0^2*h^6 + 4*RxL*D*Rx0^2*h^5 + 8*RxL*R0y*D*Rx0*h^5 - 8*RxL*D^2*Rx0*h^4*nu^2 + 16*RxL*D^2*Rx0*h^4 - 14*R0y*D^2*Rx0*h^2*nu^2 + 28*R0y*D^2*Rx0*h^2*nu + 24*R0y*D^2*Rx0*h^2 - 20*D^3*Rx0*h*nu^2 + 40*D^3*Rx0*h*nu + 48*D^3*Rx0*h + 8*RxL*R0y*D^2*h^4 - 16*RxL*D^3*h^3*nu^2 + 16*RxL*D^3*h^3 - 28*R0y*D^3*h*nu^2 + 56*R0y*D^3*h*nu + 48*R0y*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20);;
        % B1uN = ((2*K0y*RxL*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*RxL*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 12*RxL*D^2*R0y*h^2*nu^2 + 24*RxL*D^2*R0y*h^2*nu + 24*RxL*D^2*R0y*h^2 + 16*D^3*R0y*h*nu^3 - 56*D^3*R0y*h*nu^2 + 48*D^3*R0y*h + 8*K0y*RxL*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 24*RxL*D^3*h*nu^2 + 48*RxL*D^3*h*nu + 48*RxL*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)) - (8*(- 8*D^2*nu^2 + 4*RxL*h*D*nu + 8*D^2 + 4*RxL*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) - (4*D*nu)/(2*D + R0y*h) - (4*D*nu)/(2*D + RxL*h) - (2*(8*D^2*nu + 2*D*R0y*h*nu + 2*D*RxL*h*nu))/((2*D + R0y*h)*(2*D + RxL*h)) - (8*(- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) + (2*KxL*R0y*RxL^2*h^6 + 4*KxL*D*RxL^2*h^5 + 8*KxL*R0y*D*RxL*h^5 - 8*KxL*D^2*RxL*h^4*nu^2 + 16*KxL*D^2*RxL*h^4 - 12*R0y*D^2*RxL*h^2*nu^2 + 24*R0y*D^2*RxL*h^2*nu + 24*R0y*D^2*RxL*h^2 + 16*D^3*RxL*h*nu^3 - 56*D^3*RxL*h*nu^2 + 48*D^3*RxL*h + 8*KxL*R0y*D^2*h^4 - 16*KxL*D^3*h^3*nu^2 + 16*KxL*D^3*h^3 - 24*R0y*D^3*h*nu^2 + 48*R0y*D^3*h*nu + 48*R0y*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + RxL*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)) + 20);;

        % B2u0 = ((2*Kx0*Rx0*h^4 + 4*Kx0*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + Rx0*h)) - (8*D*nu)/(2*D + Rx0*h) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) + 20);;
        % B2u1 = (20 - (2*D - Rx0*h)/(2*D + Rx0*h));
        % B2u2 = 20;
        % B2uN1 = (20 - (2*D - RxL*h)/(2*D + Rx0*h));
        % B2uN = ((2*KxL*RxL*h^4 + 4*KxL*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + RxL*h)) - (8*D*nu)/(2*D + RxL*h) - (8*(4*D + 4*D*nu))/(2*D + RxL*h) + 20);;

        % BN1u0  = ((4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) - (4*D*nu)/(2*D + Rx0*h) + (2*RLy*R0y*Rx0^2*h^6 + 4*RLy*D*Rx0^2*h^5 + 8*RLy*R0y*D*Rx0*h^5 - 8*RLy*D^2*Rx0*h^4*nu^2 + 16*RLy*D^2*Rx0*h^4 - 14*R0y*D^2*Rx0*h^2*nu^2 + 28*R0y*D^2*Rx0*h^2*nu + 24*R0y*D^2*Rx0*h^2 - 20*D^3*Rx0*h*nu^2 + 40*D^3*Rx0*h*nu + 48*D^3*Rx0*h + 8*RLy*R0y*D^2*h^4 - 16*RLy*D^3*h^3*nu^2 + 16*RLy*D^3*h^3 - 28*R0y*D^3*h*nu^2 + 56*R0y*D^3*h*nu + 48*R0y*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20);;
        % BN1u1  = (20 - (2*D - Rx0*h)/(2*D + Rx0*h) - (2*D - R0y*h)/(2*D + R0y*h));
        % BN1u2  = (20 - (2*D - R0y*h)/(2*D + R0y*h));
        % BN1uN1 = ((4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) - (4*D*nu)/(2*D + Rx0*h) + (2*RxL*R0y*Rx0^2*h^6 + 4*RxL*D*Rx0^2*h^5 + 8*RxL*R0y*D*Rx0*h^5 - 8*RxL*D^2*Rx0*h^4*nu^2 + 16*RxL*D^2*Rx0*h^4 - 14*R0y*D^2*Rx0*h^2*nu^2 + 28*R0y*D^2*Rx0*h^2*nu + 24*R0y*D^2*Rx0*h^2 - 20*D^3*Rx0*h*nu^2 + 40*D^3*Rx0*h*nu + 48*D^3*Rx0*h + 8*RxL*R0y*D^2*h^4 - 16*RxL*D^3*h^3*nu^2 + 16*RxL*D^3*h^3 - 28*R0y*D^3*h*nu^2 + 56*R0y*D^3*h*nu + 48*R0y*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20);;
        % BN1uN = ((2*K0y*RxL*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*RxL*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 12*RxL*D^2*R0y*h^2*nu^2 + 24*RxL*D^2*R0y*h^2*nu + 24*RxL*D^2*R0y*h^2 + 16*D^3*R0y*h*nu^3 - 56*D^3*R0y*h*nu^2 + 48*D^3*R0y*h + 8*K0y*RxL*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 24*RxL*D^3*h*nu^2 + 48*RxL*D^3*h*nu + 48*RxL*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)) - (8*(- 8*D^2*nu^2 + 4*RxL*h*D*nu + 8*D^2 + 4*RxL*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) - (4*D*nu)/(2*D + R0y*h) - (4*D*nu)/(2*D + RxL*h) - (2*(8*D^2*nu + 2*D*R0y*h*nu + 2*D*RxL*h*nu))/((2*D + R0y*h)*(2*D + RxL*h)) - (8*(- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2) + (2*KxL*R0y*RxL^2*h^6 + 4*KxL*D*RxL^2*h^5 + 8*KxL*R0y*D*RxL*h^5 - 8*KxL*D^2*RxL*h^4*nu^2 + 16*KxL*D^2*RxL*h^4 - 12*R0y*D^2*RxL*h^2*nu^2 + 24*R0y*D^2*RxL*h^2*nu + 24*R0y*D^2*RxL*h^2 + 16*D^3*RxL*h*nu^3 - 56*D^3*RxL*h*nu^2 + 48*D^3*RxL*h + 8*KxL*R0y*D^2*h^4 - 16*KxL*D^3*h^3*nu^2 + 16*KxL*D^3*h^3 - 24*R0y*D^3*h*nu^2 + 48*R0y*D^3*h*nu + 48*R0y*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + RxL*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*RxL*h + R0y*RxL*h^2)) + 20);;

        % BNu0 =  ((2*KLy*Rx0*RLy^2*h^6 + 4*KLy*D*RLy^2*h^5 + 8*KLy*Rx0*D*RLy*h^5 - 8*KLy*D^2*RLy*h^4*nu^2 + 16*KLy*D^2*RLy*h^4 - 12*Rx0*D^2*RLy*h^2*nu^2 + 24*Rx0*D^2*RLy*h^2*nu + 24*Rx0*D^2*RLy*h^2 + 16*D^3*RLy*h*nu^3 - 56*D^3*RLy*h*nu^2 + 48*D^3*RLy*h + 8*KLy*Rx0*D^2*h^4 - 16*KLy*D^3*h^3*nu^2 + 16*KLy*D^3*h^3 - 24*Rx0*D^3*h*nu^2 + 48*Rx0*D^3*h*nu + 48*Rx0*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + RLy*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2)) - (8*(- 8*D^2*nu^2 + 4*Rx0*h*D*nu + 8*D^2 + 4*Rx0*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2) - (4*D*nu)/(2*D + RLy*h) - (4*D*nu)/(2*D + Rx0*h) - (2*(8*D^2*nu + 2*D*RLy*h*nu + 2*D*Rx0*h*nu))/((2*D + RLy*h)*(2*D + Rx0*h)) - (8*(- 8*D^2*nu^2 + 4*RLy*h*D*nu + 8*D^2 + 4*RLy*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2) + (2*Kx0*RLy*Rx0^2*h^6 + 4*Kx0*D*Rx0^2*h^5 + 8*Kx0*RLy*D*Rx0*h^5 - 8*Kx0*D^2*Rx0*h^4*nu^2 + 16*Kx0*D^2*Rx0*h^4 - 12*RLy*D^2*Rx0*h^2*nu^2 + 24*RLy*D^2*Rx0*h^2*nu + 24*RLy*D^2*Rx0*h^2 + 16*D^3*Rx0*h*nu^3 - 56*D^3*Rx0*h*nu^2 + 48*D^3*Rx0*h + 8*Kx0*RLy*D^2*h^4 - 16*Kx0*D^3*h^3*nu^2 + 16*Kx0*D^3*h^3 - 24*RLy*D^3*h*nu^2 + 48*RLy*D^3*h*nu + 48*RLy*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2)) + 20);;
        % BNu1 =  ((4*D^2*nu^2 - 4*D^2 - 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + RLy*h) - (4*D*nu)/(2*D + RLy*h) + (2*KLy*Rx0*RLy^2*h^6 + 4*KLy*D*RLy^2*h^5 + 8*KLy*Rx0*D*RLy*h^5 - 8*KLy*D^2*RLy*h^4*nu^2 + 16*KLy*D^2*RLy*h^4 - 14*Rx0*D^2*RLy*h^2*nu^2 + 28*Rx0*D^2*RLy*h^2*nu + 24*Rx0*D^2*RLy*h^2 - 20*D^3*RLy*h*nu^2 + 40*D^3*RLy*h*nu + 48*D^3*RLy*h + 8*KLy*Rx0*D^2*h^4 - 16*KLy*D^3*h^3*nu^2 + 16*KLy*D^3*h^3 - 28*Rx0*D^3*h*nu^2 + 56*Rx0*D^3*h*nu + 48*Rx0*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + RLy*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2)) - (8*D*Rx0*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*Rx0*h + RLy*Rx0*h^2) + 20);;
        % BNu2 =  ((2*KLy*RLy*h^4 + 4*KLy*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + RLy*h)) - (8*D*nu)/(2*D + RLy*h) - (8*(4*D + 4*D*nu))/(2*D + RLy*h) + 20);;
        % BNuN1 =  ((4*D^2*nu^2 - 4*D^2 - 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2) - (8*(4*D + 4*D*nu))/(2*D + RLy*h) - (4*D*nu)/(2*D + RLy*h) + (2*KLy*RxL*RLy^2*h^6 + 4*KLy*D*RLy^2*h^5 + 8*KLy*RxL*D*RLy*h^5 - 8*KLy*D^2*RLy*h^4*nu^2 + 16*KLy*D^2*RLy*h^4 - 14*RxL*D^2*RLy*h^2*nu^2 + 28*RxL*D^2*RLy*h^2*nu + 24*RxL*D^2*RLy*h^2 - 20*D^3*RLy*h*nu^2 + 40*D^3*RLy*h*nu + 48*D^3*RLy*h + 8*KLy*RxL*D^2*h^4 - 16*KLy*D^3*h^3*nu^2 + 16*KLy*D^3*h^3 - 28*RxL*D^3*h*nu^2 + 56*RxL*D^3*h*nu + 48*RxL*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + RLy*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2)) - (8*D*RxL*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2) + 20);;
        % BNuN =  ((2*KLy*RxL*RLy^2*h^6 + 4*KLy*D*RLy^2*h^5 + 8*KLy*RxL*D*RLy*h^5 - 8*KLy*D^2*RLy*h^4*nu^2 + 16*KLy*D^2*RLy*h^4 - 12*RxL*D^2*RLy*h^2*nu^2 + 24*RxL*D^2*RLy*h^2*nu + 24*RxL*D^2*RLy*h^2 + 16*D^3*RLy*h*nu^3 - 56*D^3*RLy*h*nu^2 + 48*D^3*RLy*h + 8*KLy*RxL*D^2*h^4 - 16*KLy*D^3*h^3*nu^2 + 16*KLy*D^3*h^3 - 24*RxL*D^3*h*nu^2 + 48*RxL*D^3*h*nu + 48*RxL*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + RLy*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2)) - (8*(- 8*D^2*nu^2 + 4*RxL*h*D*nu + 8*D^2 + 4*RxL*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2) - (4*D*nu)/(2*D + RLy*h) - (4*D*nu)/(2*D + RxL*h) - (2*(8*D^2*nu + 2*D*RLy*h*nu + 2*D*RxL*h*nu))/((2*D + RLy*h)*(2*D + RxL*h)) - (8*(- 8*D^2*nu^2 + 4*RLy*h*D*nu + 8*D^2 + 4*RLy*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2) + (2*KxL*RLy*RxL^2*h^6 + 4*KxL*D*RxL^2*h^5 + 8*KxL*RLy*D*RxL*h^5 - 8*KxL*D^2*RxL*h^4*nu^2 + 16*KxL*D^2*RxL*h^4 - 12*RLy*D^2*RxL*h^2*nu^2 + 24*RLy*D^2*RxL*h^2*nu + 24*RLy*D^2*RxL*h^2 + 16*D^3*RxL*h*nu^3 - 56*D^3*RxL*h*nu^2 + 48*D^3*RxL*h + 8*KxL*RLy*D^2*h^4 - 16*KxL*D^3*h^3*nu^2 + 16*KxL*D^3*h^3 - 24*RLy*D^3*h*nu^2 + 48*RLy*D^3*h*nu + 48*RLy*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + RxL*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*RLy*h + 2*D*RxL*h + RLy*RxL*h^2)) + 20);;
        % else
        [B0u0,~,~,~,~,~]          = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,0);
        [B0u1,~,~,~,~,~,~,~]      = D01_coeffs(K0y,R0y,Rx0,h,D,nu,0);
        [B0u2,~,~,~,~,~,~,~,~]    = D02_coeffs(K0y,R0y,h,D,nu,0);
        [B0uN1,~,~,~,~,~,~,~]  = D01_coeffs(K0y,R0y,RxL,h,D,nu,0);
        [B0uN,~,~,~,~,~]          = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,0);

        [B1u0,~,~,~,~,~,~,~]           = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,0);
        [B1u1,~,~,~,~,~,~,~,~,~,~]     = D11_coeffs(R0y,Rx0,h,D,nu,0);
        [B1u2,~,~,~,~,~,~,~,~,~,~,~]   = D12_coeffs(R0y,h,D,nu,0);
        [B1uN1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu,0);
        [B1uN,~,~,~,~,~,~,~]           = D10_coeffs(R0y,KxL,RxL,h,D,nu,0);

        [B2u0,~,~,~,~,~,~,~,~]           = D20_coeffs(Kx0,Rx0,h,D,nu,0);
        [B2u1,~,~,~,~,~,~,~,~,~,~,~]     = D21_coeffs(Rx0,h,D,nu,0);
        [~,~,~,~,~,~,B2u2,~,~,~,~,~,~]   = D22_coeffs;
        [B2uN1,~,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu,0);
        [B2uN,~,~,~,~,~,~,~,~]           = D20_coeffs(KxL,RxL,h,D,nu,0);

        [BN1u0,~,~,~,~,~,~,~]           = D10_coeffs(RLy,Kx0,Rx0,h,D,nu,0);
        [BN1u1,~,~,~,~,~,~,~,~,~,~]     = D11_coeffs(RLy,Rx0,h,D,nu,0);
        [BN1u2,~,~,~,~,~,~,~,~,~,~,~]   = D12_coeffs(RLy,h,D,nu,0);
        [BN1uN1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu,0);
        [BN1uN,~,~,~,~,~,~,~]           = D10_coeffs(RLy,KxL,RxL,h,D,nu,0);

        [BNu0,~,~,~,~,~]         = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,0);
        [BNu1,~,~,~,~,~,~,~]     = D01_coeffs(KLy,RLy,Rx0,h,D,nu,0);
        [BNu2,~,~,~,~,~,~,~,~]   = D02_coeffs(KLy,RLy,h,D,nu,0);
        [BNuN1,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu,0);
        [BNuN,~,~,~,~,~]         = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu,0);

        d0  = [B0u0, B0u1, B0u2, B0uN1, B0uN];
        d1  = [B1u0, B1u1, B1u2, B1uN1, B1uN];
        dc  = [B2u0, B2u1, B2u2, B2uN1, B2uN];
        dM1 = [BN1u0,BN1u1,BN1u2,BN1uN1,BN1uN];
        dM  = [BNu0, BNu1, BNu2, BNuN1, BNuN];

    case '1'
        [~,~,~,D00u01,~,~]            = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,3);
        [~,~,~,~,D01u02,~,~,~]        = D01_coeffs(K0y,R0y,Rx0,h,D,nu,4);
        [~,~,~,~,D02u03,~,~,~,~]      = D02_coeffs(K0y,R0y,h,D,nu,4);
        [~,~,~,D0Nm1u0N,~,~,~,~]      = D01_coeffs(K0y,R0y,RxL,h,D,nu,3);
        [~,~,~,~,~,~]                 = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu);

        [~,~,~,~, D10u11,~,~,~]  = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,4);
        [~,D11u12,~,~,~,~,~,~,~,~,~]   = D11_coeffs(R0y,Rx0,h,D,nu,1);
        [~,D12u13,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu,1);
        [~,~,~,D1Nm1u1N,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu,3);
        [~,~,~,~,~,~,~,~]       = D10_coeffs(R0y,KxL,RxL,h,D,nu);

        [~,D20u21,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu,1);
        [~,D21u22,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu,1);
        [~,~,~,~,~,~,~,~,~,~,D22u23,~,~] = D22_coeffs;
        [~,~,~,D2Nm1u2N,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu,3);
        [~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu);

        [~,~,~,~, DN10u11,~,~,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu,4);
        [~,DN11u12,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu,1);
        [~,DN12u13,~,~,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu,1);
        [~,~,~,DN1Nm1u1N,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu,3);
        [~,~,~,~,~,~,~,~] = D10_coeffs(RLy,KxL,RxL,h,D,nu);

        [~,~,~,DN0u01,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,3);
        [~,~,~,~,DN1u02,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu,4);
        [~,~,~,~,DN2u03,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu,4);
        [~,~,~,DNNm1u0N,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu,3);
        [~,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu);

        d0  = [D00u01;D01u02;D02u03;D0Nm1u0N];
        d1  = [D10u11;D11u12;D12u13;D1Nm1u1N];
        dc  = [D20u21;D21u22;D22u23;D2Nm1u2N];
        dM1 = [DN10u11;DN11u12;DN12u13;DN1Nm1u1N];
        dM  = [DN0u01;DN1u02;DN2u03;DNNm1u0N];

    case '2'

        [~,~,~,~,D00u02,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,4);
        [~,~,~,~,~,D01u03,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu,5);
        [~,~,~,~,~,D02u04,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu,5);

        [~,~,~,~,~, D10u12,~,~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,5);
        [~,~,D11u13,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu,2);
        [~,~,D12u14,~,~,~,~,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu,2);

        [~,~,D20u22,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu,2);
        [~,~,D21u23,~,~,~,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu,2);
        [~,~,~,~,~,~,~,~,~,~,~,~,D22u24] = D22_coeffs;

        [~,~,~,~,~, DN10u12,~,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu,5);
        [~,~,DN11u13,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu,2);
        [~,~,DN12u14,~,~,~,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu,2);

        [~,~,~,~,DN0u02,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu,4);
        [~,~,~,~,~,DN1u03,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu,5);
        [~,~,~,~,~,DN2u04,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu,5);

        d0 = [D00u02;D01u03;D02u04];
        d1 = [D10u12;D11u13;D12u14];
        dc = [D20u22;D21u23;D22u24];
        dM1 = [DN10u12;DN11u13;DN12u14];
        dM = [DN0u02;DN1u03;DN2u04];

    case 'ny-1'

        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,D01u10] = D01_coeffs(K0y,R0y,Rx0,h,D,nu,7);
        [~,~,~,~,~,~,~,~,D02u11] = D02_coeffs(K0y,R0y,h,D,nu,8);
        [~,~,~,~,~,~,D0Nm1u1Nm2,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu,6);
        [~,~,~,~,~,D0Nu1Nm1] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,5);

        [~,~,~,~,~,~,~,~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,D11u20,~,~] = D11_coeffs(R0y,Rx0,h,D,nu,8);
        [~,~,~,~,~,~,~,~,~,D12u21,~,~] = D12_coeffs(R0y,h,D,nu,8);
        [~,~,~,~,~,~,~,D1Nm1u2Nm2,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu);
        [~,~,~,~,~,~, D1Nu2Nm1,~] = D10_coeffs(R0y,KxL,RxL,h,D,nu);

        [~,~,~,~,~,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,~,D21u30,~,~] = D21_coeffs(Rx0,h,D,nu);
        [~,~,~,D22u31,~,~,~,~,~,~,~,~,~] = D22_coeffs;
        [~,~,~,~,~,~,~,~,D2Nm1u3Nm2,~,~,~] = D21_coeffs(RxL,h,D,nu);
        [~,~,~,~,~,~,~,D2Nu3Nm1,~] = D20_coeffs(KxL,RxL,h,D,nu);

        [~,~,~,~,~,~,~,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,DN11u00,~] = D11_coeffs(RLy,Rx0,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,DN12u01,~] = D12_coeffs(RLy,h,D,nu) ;
        [~,~,~,~,~,~,~,~,~,~,DN1Nm1u0Nm2] = D11_coeffs(RLy,RxL,h,D,nu) ;
        [~,~,~,~,~,~,~,DN1Nu0Nm1] = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;

        d0 = [D01u10;D02u11;D0Nm1u1Nm2;D0Nu1Nm1];
        d1 = [D11u20;D12u21;D1Nm1u2Nm2;D1Nu2Nm1];
        dc = [D21u30;D22u31;D2Nm1u3Nm2;D2Nu3Nm1];
        dM1 =[DN11u00;DN12u01;DN1Nm1u0Nm2;DN1Nu0Nm1];
        dM = [NaN; NaN; NaN; NaN;];

    case 'ny'

        [~,D00u10,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,1);
        [~,D01u11,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu,1);
        [~,D02u12,~,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu,1);
        [~,D0Nm1u1Nm1,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu,1);
        [~,D0Nu1N,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,1);

        [~, D10u20,~,~,~,~,~,~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,1);
        [~,~,~,~,~,D11u21,~,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu,5);
        [~,~,~,~,~,~,D12u22,~,~,~,~,~] = D12_coeffs(R0y,h,D,nu,6);
        [~,~,~,~,~,D1Nm1u2Nm1,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu,5);
        [~, D1Nu2N,~,~,~,~,~,~] = D10_coeffs(R0y,KxL,RxL,h,D,nu,1);

        [~,~,~,~,D20u30,~,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu,4);
        [~,~,~,~,~,D21u31,~,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu,5);
        [~,~,~,~,~,~,~,D22u32,~,~,~,~,~] = D22_coeffs;
        [~,~,~,~,~,D2Nm1u3Nm1,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu,5);
        [~,~,~,~,D2Nu3N,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu,4);

        [~,~,~, DN10u00,~,~,~,~] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu,3);
        [~,~,~,~,DN11u01,~,~,~,~,~,~] = D11_coeffs(RLy,Rx0,h,D,nu,4);
        [~,~,~,~,~,DN12u02,~,~,~,~,~,~] = D12_coeffs(RLy,h,D,nu,5);
        [~,~,~,~,DN1Nm1u0Nm1,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu,4);
        [~,~,~, DN1Nu0N,~,~,~,~] = D10_coeffs(RLy,KxL,RxL,h,D,nu,3);

        d0 = [D00u10;D01u11;D02u12;D0Nm1u1Nm1;D0Nu1N];
        d1 = [D10u20;D11u21;D12u22;D1Nm1u2Nm1;D1Nu2N];
        dc = [D20u30;D21u31;D22u32;D2Nm1u3Nm1;D2Nu3N];
        dM1 = [DN10u00;DN11u01;DN12u02;DN1Nm1u0Nm1;DN1Nu0N];
        dM = [NaN; NaN; NaN; NaN;NaN];

    case 'ny+1'

        [~,~,~,~,~,D00u11] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,5);
        [~,~,~,~,~,~,D01u12,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu,6);
        [~,~,~,~,~,~,~,D02u13,~] = D02_coeffs(K0y,R0y,h,D,nu,7);
        [~,~,~,~,~,~,~,D0Nm1u1N] = D01_coeffs(K0y,R0y,RxL,h,D,nu,7);
        [~,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu);

        [~,~,~,~,~,~, D10u21,~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,6);
        [~,~,~,~,~,~,~,D11u22,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu,7);
        [~,~,~,~,~,~,~,~,D12u23,~,~,~] = D12_coeffs(R0y,h,D,nu,8);
        [~,~,~,~,~,~,~,~,D1Nm1u2N,~,~] = D11_coeffs(R0y,RxL,h,D,nu,8);
        [~,~,~,~,~,~,~,~] = D10_coeffs(R0y,KxL,RxL,h,D,nu);

        [~,~,~,~,~,~,~,D20u31,~] = D20_coeffs(Kx0,Rx0,h,D,nu,7);
        [~,~,~,~,~,~,~,~,D21u32,~,~,~] = D21_coeffs(Rx0,h,D,nu);
        [~,~,~,~,~,~,~,~,~,~,~,D22u33,~] = D22_coeffs;
        [~,~,~,~,~,~,~,~,~,D2Nm1u3N,~,~] = D21_coeffs(RxL,h,D,nu);
        [~,~,~,~,~,~,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu);

        [~,~,~,~,~,~,~, DN10u01] = D10_coeffs(RLy,Kx0,Rx0,h,D,nu,7);
        [~,~,~,~,~,~,~,~,~,~,DN11u02] = D11_coeffs(RLy,Rx0,h,D,nu,8);
        [~,~,~,~,~,~,~,~,~,~,~,DN12u03] = D12_coeffs(RLy,h,D,nu);
        [~,~,~,~,~,~,~,~,~,DN1Nm1u0N,~] = D11_coeffs(RLy,RxL,h,D,nu,8);
        [~,~,~,~,~,~,~,~] = D10_coeffs(RLy,KxL,RxL,h,D,nu);

        d0 =[D00u11;D01u12;D02u13;D0Nm1u1N];
        d1 =[D10u21;D11u22;D12u23;D1Nm1u2N];
        dc =[D20u31;D21u32;D22u33;D2Nm1u3N];
        dM1 =[DN10u01;DN11u02;DN12u03;DN1Nm1u0N];
        dM = [NaN; NaN; NaN; NaN;];

    case '2ny'

        [~,~,D00u20,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu,2);
        [~,~,D01u21,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu,2);
        [~,~,D02u22,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu,2);
        [~,~,D0Nm1u2Nm1,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu,2);
        [~,~,D0Nu2N,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu,2);

        [~,~, D10u30,~,~,~,~,~] = D10_coeffs(R0y,Kx0,Rx0,h,D,nu,2);
        [~,~,~,~,~,~,D11u31,~,~,~,~] = D11_coeffs(R0y,Rx0,h,D,nu,6);
        [~,~,~,~,~,~,~,D12u32,~,~,~,~] = D12_coeffs(R0y,h,D,nu,7);
        [~,~,~,~,~,~,D1Nm1u3Nm1,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu,6);
        [~,~, D1Nu3N,~,~,~,~,~] = D10_coeffs(R0y,KxL,RxL,h,D,nu,2);

        [~,~,~,~,~,D20u40,~,~,~] = D20_coeffs(Kx0,Rx0,h,D,nu,5);
        [~,~,~,~,~,~,D21u41,~,~,~,~,~] = D21_coeffs(Rx0,h,D,nu,6);
        [~,~,~,~,~,~,~,~,D22u42,~,~,~,~] = D22_coeffs;
        [~,~,~,~,~,~,D2Nm1u4Nm1,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu,6);
        [~,~,~,~,~,D2Nu4N,~,~,~] = D20_coeffs(KxL,RxL,h,D,nu,5);

        d0 = [D00u20;D01u21;D02u22;D0Nm1u2Nm1;D0Nu2N];
        d1 = [D10u30;D11u31;D12u32;D1Nm1u3Nm1;D1Nu3N];
        dc = [D20u40;D21u41;D22u42;D2Nm1u4Nm1;D2Nu4N];
        dM1 = [NaN; NaN; NaN; NaN;NaN];
        dM = [NaN; NaN; NaN; NaN;NaN];

end


end

