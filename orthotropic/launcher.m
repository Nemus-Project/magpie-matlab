clear all 
close all
clc

%-------------------------------------------------------------------------
%               MAGPIE launcher
%                Dr M Ducceschi
%             University of Bologna
%                  11 Sep 2024
%-------------------------------------------------------------------------

%-- Input constants

%-- NB To recover the isotropic case, for input E, nu, USE:
% Ex  = E 
% Ey  = E; 
% Gxy = E/(2*(1+nu)); 
% nux = nu ;

rho = 390 ;
Ex  = 10.9e9 ;
Ey  = 0.64e9 ;
Gxy = 0.58e9 ;
nux = 0.39 ;

Lx  = 0.6 ;
Ly  = 0.4 ;
Lz  = 1e-3 ;

hx  = Lx/100 ;
hy  = Ly/100 ;

%-- elastic constants around the edges
K0y     = 1e15 ;
R0y     = 1e15 ;
Kx0     = 0e15 ;
Rx0     = 0e15 ;
KLy     = 0e15 ;
RLy     = 0e15;
KxL     = 0e15 ;
RxL     = 0e15 ;

Nmodes  = 9 ;
printfreqs = 1 ;

orthoPlateEigs(rho,Ex,Ey,Gxy,nux,Lx,Ly,Lz,hx,hy,K0y,R0y,Kx0,Rx0,KLy,RLy,KxL,RxL,Nmodes,printfreqs) ;
