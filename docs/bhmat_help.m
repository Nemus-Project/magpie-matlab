%% bhmat 
%
% Generate the 2D Biharmonic operator with elastic boundary constants.
%
%% Syntax
%
%   biharm = bhmat(BCs,Nxy,h,Lz,E,nu)
%
%% Description
%
% |biharm = bhmat(BCs,Nxy,h,Lz,E,nu)| returns a |Nxy|-by-|Nxy|  discrete biharmonic 
% operator matrix |biharm| with elastic constants |BCs|, grid spacing |h|, plate thickness
% |Lz|, youngs Modulus |E| and Poisson's number  |nu|
%
%% Example
%
%   %% physical and elastic parameters
%   Lz  = 0.8e-3;
%   E   = 1.01e+11;
%   nu  = 0.3;
%   Nx  = 100;
%   Ny  = 100;
%   h   = 1e-2;
%   BCs = 1e15 * ones(4,2) %-- Clamped condition on all sides
%   
%   biharm = bhmat(BCs,[Nx Ny], h, Lz, E, nu);
%
%% Input Arguments
% * |BCs|  : 4-by-2 matrix of elastic boundary constants, in the order y0, x0, yL, xL. Colun 1 is a flexural constraint and column 2 a rotational constraint
% * |Nxy|  : 2 element vector with number of gridpoints in each dimension |[Nx Ny]|
% * |h|    : Grid spacing in metres
% * |Lz|   : Plate thickness in metres
% * |E|    : Young's mod [Pa]
% * |nu|   : Poisson's number
%
%% See Also
% |magpie|

