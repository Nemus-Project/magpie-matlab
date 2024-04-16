%% BHMAT Generate 2D Biharmonic Operator
%   biHarm = BHMAT(BCs,Nxy,h,Lz,E,nu) % compute a biharmonic for a  2D Thin Plate
%
%% Arguments:
%
%       BCS          %-- density [kg/m^3]
%
% Boundary condition columns represent a different elastic stiffness
% contant for each side of the plate
%
%  - Column 1 flexural constant
%  - Column 2 rotational constant
%
%       BCs = [K0y, R0y;
%              Kx0, Rx0;
%              KLy, RLy;
%              KxL, RxL];
%
%       Nxy          %-- Number of grid points [Nx Ny]
%       h            %-- Grid spacing
%       Lz           %-- plate thickness
%       E            %-- Young's mod [Pa]
%       nu           %-- poisson's ratio
%
%% Returns
%
% Biharmonic operator with size [Nx Ny]
function biharm = bhmat(BCs,Nxy,h,Lz,E,nu)

%% Argument Validation
validateattributes(BCs,      {'double'}, {'size', [4,2]});
validateattributes(Nxy,      {'double'}, {'numel', 2, 'positive','integer'});
validateattributes(h,        {'double'}, {'nonempty'});
validateattributes(Lz,       {'double'}, {'nonempty'});
validateattributes(E,        {'double'}, {'nonempty'});
validateattributes(nu,       {'double'}, {'nonempty'});

%% Unpack Variables
%
Nx = Nxy(1);
Ny = Nxy(2);

%% Define Commonly used Variables
D = E * Lz^3 / 12 / (1-nu^2);

a0 = ones(Ny+1,1) ;
a1 = ones(Ny,1) ;
a2 = ones(Ny-1,1) ;

biharm = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
%% Define Diagonals
%%%% dm2Ny  

[~,~,dc,dM1,dM] = biharmdiag(BCs, h, D, nu,'-2ny');

dm2Ny0 = dc(3)*a0;
dm2Ny0([1,2,Ny,Ny+1]) = [dc(1),dc(2),dc(end-1),dc(end)];

dm2Ny1 = dM1(3)*a0;
dm2Ny1([1,2,Ny,Ny+1]) = [dM1(1),dM1(2),dM1(end-1),dM1(end)];

dm2Ny2 = dM(3)*a0;
dm2Ny2([1,2,Ny,Ny+1]) = [dM(1),dM(2),dM(end-1),dM(end)];

Dm2Ny = [repmat(dm2Ny0,Nx-3,1);dm2Ny1;dm2Ny2];

%% dmNym1 

[~,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'-ny-1');

dmNym10 = d1(2)*a1 ;
dmNym10([1,Ny-1,Ny]) = [d1(1), d1(end-1), d1(end)];
dmNym11 = dc(2)*a1 ;
dmNym11([1,Ny-1,Ny]) = [dc(1), dc(end-1), dc(end)];
dmNym1M1 = dM1(2)*a1;
dmNym1M1([1,Ny-1,Ny]) = [dM1(1), dM1(end-1), dM1(end)];
dmNym1M = dM(2)*a1;
dmNym1M([1,Ny-1,Ny]) = [dM(1), dM(end-1), dM(end)];

DmNym1 = [dmNym10;0;repmat([dmNym11;0],Nx-3,1);dmNym1M1;0;dmNym1M];

%% dmNy   

[~,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'-ny');

dmNy0 = d1(3)*a0;
dmNy0([1,2,Ny,Ny+1]) = [d1(1),d1(2),d1(end-1),d1(end)];
dmNy1 = dc(3)*a0;
dmNy1([1,2,Ny,Ny+1]) = [dc(1),dc(2),dc(end-1),dc(end)];
dmNyM1 = dM1(3)*a0;
dmNyM1([1,2,Ny,Ny+1]) = [dM1(1),dM1(2),dM1(end-1),dM1(end)];
dmNyM = dM(3)*a0;
dmNyM([1,2,Ny,Ny+1]) = [dM(1),dM(2),dM(end-1),dM(end)];

DmNy = [dmNy0;repmat(dmNy1,(Nx-3),1);dmNyM1;dmNyM];

%% dmNyp1 % pad zeros at the

[~,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'-ny+1');

dmNy10 = d1(2)*a1 ;
dmNy10([1,2,Ny]) = [d1(1), d1(2), d1(end)];
dmNy11 = dc(2)*a1 ;
dmNy11([1,2,Ny]) = [dc(1), dc(2), dc(end)];
dmNy1M1 = dM1(2)*a1;
dmNy1M1([1,2,Ny]) = [dM1(1), dM1(2), dM1(end)];
dmNy1M = dM(2)*a1;
dmNy1M([1,2,Ny]) = [dM(1), dM(2), dM(end)];

DmNy1 = [0;dmNy10;repmat([0;dmNy11],Nx-3,1);0;dmNy1M1;0;dmNy1M;0];

%% dm2   
[d0,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'-2');

dm20 = d0(1)*a2;
dm20([Ny-2,Ny-1]) = [d0(end-1),d0(end)];

dm21 = d1(1)*a2;
dm21([Ny-2,Ny-1]) = [d1(end-1),d1(end)];

dm22 = dc(1)*a2;
dm22([Ny-2,Ny-1]) = [dc(end-1),dc(end)];

dm2M1 = dM1(1)*a2;
dm2M1([Ny-2,Ny-1]) = [dM1(end-1),dM1(end)];

dm2M = dM(1)*a2;
dm2M([Ny-2,Ny-1]) = [dM(end-1),dM(end)];

Dm2 = [dm20;0;0;dm21;0;0;repmat([dm22;0;0],Nx+1-4,1);dm2M1;0;0;dm2M];

%% dm1   

[d0,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'-1');

dm10 = d0(2)*a1;
dm10([1,Ny-1,Ny]) = [d0(1),d0(end-1),d0(end)];
dm11 = d1(2)*a1;
dm11([1,Ny-1,Ny]) = [d1(1),d1(end-1),d1(end)];
dm12 = dc(2)*a1;
dm12([1,Ny-1,Ny]) = [dc(1),dc(end-1),dc(end)];
dm1M1 = dM1(2)*a1;
dm1M1([1,Ny-1,Ny]) = [dM1(1),dM1(end-1),dM1(end)];
dm1M = dM(2)*a1;
dm1M([1,Ny-1,Ny]) = [dM(1),dM(end-1),dM(end)];

Dm1 = [dm10;0;dm11;0;repmat([dm12;0],Nx+1-4,1);dm1M1;0;dm1M];

%% d00
[d0,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'0');

d00 = d0(3)*a0;
d00([1,2,Ny,Ny+1]) = [d0(1),d0(2),d0(end-1),d0(end)];
d01 = d1(3)*a0;
d01([1,2,Ny,Ny+1]) = [d1(1),d1(2),d1(end-1),d1(end)];
d02 = dc(3)*a0;
d02([1,2,Ny,Ny+1]) = [dc(1),dc(2),dc(end-1),dc(end)];
d0M1 = dM1(3)*a0;
d0M1([1,2,Ny,Ny+1]) = [dM1(1),dM1(2),dM1(end-1),dM1(end)];
d0M = dM(3)*a0;
d0M([1,2,Ny,Ny+1]) = [dM(1),dM(2),dM(end-1),dM(end)];

D0 = [d00;d01;repmat(d02,(Nx+1-4),1);d0M1;d0M];

%% dp1   

[d0,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'1');

d10 = d0(3)*a1;
d10([1,2,Ny]) = [d0(1),d0(2),d0(end)];
d11 = d1(3)*a1;
d11([1,2,Ny]) = [d1(1),d1(2),d1(end)];
d12 = dc(3)*a1;
d12([1,2,Ny]) = [dc(1),dc(2),dc(end)];
d1M1 = dM1(3)*a1;
d1M1([1,2,Ny]) = [dM1(1),dM1(2),dM1(end)];
d1M = dM(3)*a1;
d1M([1,2,Ny]) = [dM(1),dM(2),dM(end)];


D1 = [d10;0;d11;0;repmat([d12;0],Nx+1-4,1);d1M1;0;d1M];

%% dp2   

[d0,d1,dc,dM1,dM] =biharmdiag(BCs, h, D, nu,'2');

d20 = d0(3)*a2;
d20([1,2,Ny-1]) = [d0(1),d0(2),d0(end)];

d21 = d1(3)*a2;
d21([1,2,Ny-1]) = [d1(1),d1(2),d1(end)];

d22 = dc(3)*a2;
d22([1,2,Ny-1]) = [dc(1),dc(2),dc(end)];

d2M1 = dM1(3)*a2;
d2M1([1,2,Ny-1]) = [dM1(1),dM1(2),dM1(end)];

d2M = dM(3)*a2;
d2M([1,2,Ny-1]) = [dM(1),dM(2),dM(end)];

D2 = [d20;0;0;d21;0;0;repmat([d22;0;0],Nx+1-4,1);d2M1;0;0;d2M];

%% dpNym1 

[d0,d1,dc,dM1,~] =biharmdiag(BCs, h, D, nu,'ny-1');

dpNym10 = d0(2)*a1 ;
dpNym10([1,Ny-1,Ny]) = [d0(1), d0(end-1), d0(end)];
dpNym11 = d1(2)*a1 ;
dpNym11([1,Ny-1,Ny]) = [d1(1), d1(end-1), d1(end)];
dpNym12 = dc(2)*a1;
dpNym12([1,Ny-1,Ny]) = [dc(1), dc(end-1), dc(end)];
dpNym1M = dM1(2)*a1;
dpNym1M([1,Ny-1,Ny]) = [dM1(1), dM1(end-1), dM1(end)];

DNym1 = [0;dpNym10;0;dpNym11;repmat([0;dpNym12],Nx-3,1);0;dpNym1M;0];

%% dpNy   
[d0,d1,dc,dM1,~] =biharmdiag(BCs, h, D, nu,'ny');

dpNy0 = d0(3)*a0 ;
dpNy0([1,2,Ny,Ny+1]) = [d0(1), d0(2), d0(end-1),d0(end)];
dpNy1 = d1(3)*a0 ;
dpNy1([1,2,Ny,Ny+1]) = [d1(1), d1(2), d1(end-1),d1(end)];
dpNy2 = dc(3)*a0;
dpNy2([1,2,Ny,Ny+1]) = [dc(1), dc(2), dc(end-1),dc(end)];
dpNyM = dM1(3)*a0;
dpNyM([1,2,Ny,Ny+1]) = [dM1(1), dM1(2), dM1(end-1),dM1(end)];


DNy = [dpNy0;dpNy1;repmat(dpNy2,(Nx-3),1);dpNyM];

%% dpNyp1 
[d0,d1,dc,dM1,~] =biharmdiag(BCs, h, D, nu,'ny+1');

dpNy10 = d0(3)*a1 ;
dpNy10([1,2,Ny]) = [d0(1), d0(2), d0(end)];

dpNy11 = d1(3)*a1 ;
dpNy11([1,2,Ny]) = [d1(1), d1(2), d1(end)];

dpNy12 = dc(3)*a1;
dpNy12([1,2,Ny]) = [dc(1), dc(2), dc(end)];

dpNy1M = dM1(3)*a1;
dpNy1M([1,2,Ny]) = [dM1(1), dM1(2), dM1(end)];

DNy1 = [dpNy10;0;dpNy11;0;repmat([dpNy12;0],Nx-3,1);dpNy1M];

%% dp2Ny  

[d0,d1,dc,~,~]  =biharmdiag(BCs, h, D, nu,'2ny');

dp2Ny0 = d0(3)*a0;
dp2Ny0([1,2,Ny,Ny+1]) = [d0(1),d0(2),d0(end-1),d0(end)];

dp2Ny1 = d1(3)*a0;
dp2Ny1([1,2,Ny,Ny+1]) = [d1(1),d1(2),d1(end-1),d1(end)];

dp2Ny2 = dc(3)*a0;
dp2Ny2([1,2,Ny,Ny+1]) = [dc(1),dc(2),dc(end-1),dc(end)];

D2Ny = [dp2Ny0;dp2Ny1;repmat(dp2Ny2,Nx-3,1)];

%% Zero Padding
% negative diagonals pad at the end
Dm2Ny = [Dm2Ny; zeros((Ny+1)*2,1)];

DmNym1 = [DmNym1; zeros((Ny + 2),1)];
DmNy = [DmNy; zeros((Ny + 1),1)];
DmNy1 = [DmNy1; zeros((Ny),1) ];

Dm2 = [Dm2;0;0];
Dm1 = [Dm1;0];

% While positive diagonals pad at the start

D1 = [0;D1];
D2 = [0;0;D2];

DNym1 = [zeros((Ny),1);DNym1];
DNy = [zeros((Ny+1),1);DNy];
DNy1 = [zeros((Ny+2),1);DNy1];

D2Ny = [zeros((Ny+1)*2,1);D2Ny];
%% Assembling Diagonals into the Biharmonic

BHdiags = [Dm2Ny,...
    DmNym1, DmNy, DmNy1,...
    Dm2, Dm1, D0, D1, D2,...
    DNym1, DNy, DNy1,...
    D2Ny];

dn = [-(2*(Ny+1)),...
      -(Ny+2):(-Ny),... 
       (-2:2),...
       Ny:(Ny+2),...
       2*(Ny+1)];

BH = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
biharm = ((1/h)^4) * spdiags(BHdiags, dn, BH);

end
