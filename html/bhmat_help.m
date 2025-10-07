%% BHMAT 
% Generate 2D Biharmonic Operator
%% Syntax
%   biHarm = BHMAT(BCs,Nxy,h,Lz,E,nu) % compute a biharmonic for a  2D Thin Plate
%% Arguments
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