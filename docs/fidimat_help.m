%% fidimat 
%
% Generate Finite Difference Spatial Sparcity Matrices
%
%% Syntax
%
%   fdmat = fidimat(l,ord)
%   fdmat = fidimat(l,m,ord)
%   fdmat = fidimat(l,m,ord,bctype)
%
%% Description
%
% |fdmat = fidimat(l,ord)| generates 1D stencil |fdmat| with order given by |ord|
%
% |fdmat = fidimat(l,m,ord)| generate a stencil |fdmat| with order given by |ord|
% for a 2D system of size |l|-by-|m|. Boundary condition defaults to simply supported
%
% |fdmat = fidimat(l,m,ord,bctype)| generate stencil |fdmat| with order given
% by |ord| for a 2D system of size |l|-by-|m| with specified boundary conditions bctype.
%
%% Example
%
%   Nx = 100;
%   Ny = 100;
%   XXYY = fidimat(Ny,Nx,'xxyy', 1);  % 2D second order matrix
%
%% Input Arguments
% * m       % number of total grid points X axis
% * l       % number of total grid points Y axis
% * ord     % order of the matrix (string)
% * bctype  % boundary condition type: 1: simply supported, 2: clamped
%
%   % Valid orders
%     ['x-','x+','x.','xx','xxxx',
%     'y-','y+','y.','yy','yyyy',
%     'grad','xy','xxyy','laplace','biharm','I'];
% 
%% See Also
% <./bhmat_help.html  |bhmat|>