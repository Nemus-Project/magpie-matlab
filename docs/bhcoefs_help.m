%% bhcoefs
% Biharmonic coefficients of a block order
%% Syntax
% |varargout = BHCOEFS(BCs,h,D,nu,order)|
%% Description
% |varargout = BHCOEFS(BCs,h,D,nu,order)| returns the individual coefficients
% of the isotropic biharmonic with elastic boundary constants for a biharmonic 
% operator matrix with grid spacing |h|, stiffness coefficient |D| and Poisson's 
% number |nu|. |varargout| is a variable sized list of coefficients for the given |order|.
%
%% Input Arguments
%
% * |BCs|  : 4-by-2 matrix of elastic boundary constants, in the order y0, x0, yL, xL
% * |h|    : Grid spacing
% * |D|    : Stiffness constant $D = \frac{E L_z^3}{12(1 - \nu^2)}$
% * |nu|   : Poisson's ratio
% * |order|: Index of coefficient
%
%% See Also
