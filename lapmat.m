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
function lapl = lapmat(Nxy,h,T)

%% Argument Validation
% validateattributes(BCs,      {'double'}, {'size', [4,2]});
% validateattributes(Nxy,      {'double'}, {'numel', 2, 'positive','integer'});
% validateattributes(h,        {'double'}, {'nonempty'});
% validateattributes(Lz,       {'double'}, {'nonempty'});
% validateattributes(E,        {'double'}, {'nonempty'});
% validateattributes(nu,       {'double'}, {'nonempty'});

%% Unpack Variables
%
Nx = Nxy(1);
Ny = Nxy(2);

%% Define Commonly used Variables
%D = E * Lz^3 / 12 / (1-nu^2);

a0 = ones(Ny+1,1) ;
a1 = ones(Ny,1) ;
a2 = ones(Ny-1,1) ;

lapl = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
%% Define Diagonals


%% d00

D0 = -4*ones((Nx+1)*(Ny+1),1);

%% dp1
dp1=[ones((Ny),1);0];
dp1cut=ones((Ny),1);

D1 = [repmat(dp1,Nx,1);dp1cut];
%% dm1   

Dm1 = D1;


%% dmNyp2 

DmNy2 = ones((Nx)*(Ny+1),1);



%% dpNyp2 

DNy2 = DmNy2;

%% Zero Padding
% negative diagonals pad at the end

DmNy2 = [DmNy2; zeros((Ny + 1),1)];

Dm1 = [Dm1;0];

% While positive diagonals pad at the start

D1 = [0;D1];

DNy2 = [ zeros((Ny + 1),1);DNy2];


%% Assembling Diagonals into the Biharmonic

LPdiags = [DmNy2, Dm1, D0, D1,DNy2];

dn = [-(Ny+1),(-1:1),(Ny+1)];

LP = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
lapl = ((1/h)^2)*T * spdiags(LPdiags, dn, LP);

end
