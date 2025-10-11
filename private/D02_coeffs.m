function [u02_c,u12_c,u22_c,u01_c,u03_c,u04_c,u00_c,u13_c,u11_c] = D02_coeffs(K0y,R0y,h,D,nu,index)
if (nargin < 6)
    index = -1;
end
u02_c = NaN;
u12_c = NaN;
u22_c = NaN;
u01_c = NaN;
u03_c = NaN;
u04_c = NaN;
u00_c = NaN;
u13_c = NaN;
u11_c = NaN;

if index == 0
    u02_c = ((2*K0y*R0y*h^4 + 4*K0y*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + R0y*h)) - (8*D*nu)/(2*D + R0y*h) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) + 20) ;
elseif index == 1
    u12_c = ((8*(2*D - R0y*h))/(2*D + R0y*h) + (8*D^2*nu - 24*D^2)/(D*(2*D + R0y*h)) - 8) ;
elseif index == 2
    u22_c = ((2*D^2 + R0y*h*D)/(D*(2*D + R0y*h)) + 1) ;
elseif index == 3
    u01_c = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (- 8*D^2*nu^2 + 16*D^2*nu + 8*D^2)/(D*(2*D + R0y*h)) + (16*D*nu)/(2*D + R0y*h) - 8) ;
elseif index == 4
    u03_c = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (- 8*D^2*nu^2 + 16*D^2*nu + 8*D^2)/(D*(2*D + R0y*h)) + (16*D*nu)/(2*D + R0y*h) - 8) ;;
elseif index == 5
    u04_c = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + R0y*h)) - (4*D*nu)/(2*D + R0y*h) + 1) ;
elseif index == 6
    u00_c = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + R0y*h)) - (4*D*nu)/(2*D + R0y*h) + 1);
elseif index == 7
    u13_c = (2 - (4*D^2*nu - 8*D^2)/(D*(2*D + R0y*h)) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;
elseif index == 8
    u11_c = (2 - (4*D^2*nu - 8*D^2)/(D*(2*D + R0y*h)) - (2*(2*D - R0y*h))/(2*D + R0y*h));

else
    u02_c = ((2*K0y*R0y*h^4 + 4*K0y*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + R0y*h)) - (8*D*nu)/(2*D + R0y*h) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) + 20) ;
    u12_c = ((8*(2*D - R0y*h))/(2*D + R0y*h) + (8*D^2*nu - 24*D^2)/(D*(2*D + R0y*h)) - 8) ;
    u22_c = ((2*D^2 + R0y*h*D)/(D*(2*D + R0y*h)) + 1) ;
    u01_c = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (- 8*D^2*nu^2 + 16*D^2*nu + 8*D^2)/(D*(2*D + R0y*h)) + (16*D*nu)/(2*D + R0y*h) - 8) ;
    u03_c = u01_c;
    u04_c = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + R0y*h)) - (4*D*nu)/(2*D + R0y*h) + 1) ;
    u00_c = u04_c;
    u13_c = (2 - (4*D^2*nu - 8*D^2)/(D*(2*D + R0y*h)) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;
    u11_c = u13_c;
end
end

