%% youngcalc
% Calculate Young's Modulus from measured modal frequencies 
%   E  = youngcalc(rho,ldim,h,BCs,ExpFreq,Ntrain)
%% Description
% |E  = youngcalc(rho,ldim,h,BCs,ExpFreq,Ntrain)|estimats Young's modulus |E| of a plate
% starting from a batch of experimentally measured frequencies |ExpFreq|.
%
%% Example Usage
%
%   ExpFreq = [73.2; 148; 376; 431; 559; 910] ;  %-- these are measured from a plate 
%   rho     = 8765 ;            %-- density [kg/m^3]
%   Lx      = 0.1 ;
%   Ly      = 0.08 ;
%   Lz      = 0.00081 ;
%
%   % cantilever BCs. Clamped edge along x:
%   BCs = [0,    0;
%          1e15, 1e15;
%          0,    0;
%          0,    0];
%   
%   ldim    = [Lx Ly Lz] ;
%   h       = sqrt(Lx*Ly)*0.01 ;  %-- grid spacing [m]
%   
%   E  = youngcalc(rho,ldim,h,BCs,ExpFreq,3) ;
%
%% Input Arhguments
% * |rho|     : the experimental plate density
% * |ldim|    : a 3-by-1 array containing the Lx Ly Lz dimensions
% * |h|       : the grid spacing of the FD scheme
% * |BCs|     : a 4-by-2 array containing the rigidities of the boundary supports of the experimental plate
% * |ExpFreq| : an array contaning the first Nmodes measured modal frequencies
% * |Ntrain|  : an integer. Must be <= |Nmodes|. It is the number of training modes out of the available batch
%
%% See Also
