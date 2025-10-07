%% magpie 
%
% Modal analysis of 2D Plates with Generalised Elastic Boundary Conditions
%
%% Syntax
%   Om = magpie (density, Youngs, poisson, dim, h, BCs)
%   [Om,Q,Dm] = magpie (density, Youngs, poisson, dim, h, BCs)
%   [Om,Q,Nxy] = magpie (density, Youngs, poisson, dim, h, BCs)
%   [Om,Q,Nxy,biharm] = magpie (density, Youngs, poisson, dim, h, BCs)
%   [Om,Q,Nxy,biharm,Dm] = magpie (density, Youngs, poisson, dim, h, BCs)
%
%   [Om, Q] = magpie (density, Youngs, poisson, dim, h, BCs, Number_of_modes)
%   
%   [Om,Q] = magpie (density, Youngs, poisson, dim, h, BCs, Number_of_modes, plot_type, normalisation)
%
%% Description
% |Om = magpie (density, Youngs, poisson, dim, h, BCs)| calculates all possible angular modal
% frequencies |Om| for the given plate parameters.
%%
% |Om = magpie (density, Youngs, poisson, dim, h, BCs)|
% 
% |[Om,Q,Dm] = magpie (density, Youngs, poisson, dim, h, BCs)|
% 
% |[Om,Q,Nxy] = magpie (density, Youngs, poisson, dim, h, BCs)|
% 
% |[Om,Q,Nxy,biharm] = magpie (density, Youngs, poisson, dim, h, BCs)|
% 
% |[Om,Q,Nxy,biharm,Dm] = magpie (density, Youngs, poisson, dim, h, BCs)| optional outputs of eigenvectors |Q|
% adjusted grid size |Nxy|, biharmonic |biharm| and eigenvalues |Dm|.
%
% |Om = magpie (density, Youngs, poisson, dim, h, BCs, Number_of_modes)| avaialable with all the optional outputs above, but for a given
% number of modes. Usefull for large grids where calculating all modes is not feasible with memory constraints.
%   
%% Example
%
%   Lx   = 1.10; Ly = 0.8; Lz = 5e-3;
%   ldim = [Lx Ly Lz]; % plate dimensions [x, y, z] in metres
%   
%   E   = 117e9; %-- Young's mod [Pa]
%   rho = 8765;  %-- density [kg/m^3]
%   nu  = 0.3;   %-- poisson's ratio
%   Nm  = 6;                %-- number of modes to compute
%   h   = sqrt(Lx*Ly)*1e-2; %-- Grid Spacing
%   BCs = zeros(4,2);   %-- elastic constants
%   BCs (:,1) = 1e15;   %-- Felxural restraint high (simply supported condition)
%   
%   Om = magpie(rho,E,nu,ldim, h, BCs, Nm);
%   
%% Input Arguments
%
% * |density| : plate density in kg/m^3
% * |Youngs| : Young's modulus of plate material
% * |poisson| : Poisson number of plate material
% * |dim| : 3 element array describing plate dimensions in metres of the form |[Lx Ly Lz]|
% * |h| : Grid spacing for the plate
% * |BCs| : 4-by-2 matrix of eleastic boundary constants around each edge of
%       the plate. Column 1 flexural and Column 2 the rotational constraint
%% Output
% * |Om|      : Angular modal frequencies
% * |Q|       : A matrix of column eigenvector(s)
% * |Nxy|     : 2-element array of grid points |[Nx Ny]| used for calculation. 
% * |biharm|  : Discrete biharmonic operator matrix
% * |Dm|      : the eigenvalues of the biharmonic matrix
%
%% See Also
% <./bhmat_help.html  |bhmat|> &#124; <./youngcalc_help.html |youngcalc|>