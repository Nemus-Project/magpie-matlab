%-------------------------------------------------------------------------
%             Frequency Domain Launcher 
%   Orthotropic Plates with Ribs, Static loads and Stiffeners 
%                Dr M Ducceschi
%             University of Bologna
%                  30 Sep 2024
%-------------------------------------------------------------------------

clear all
close all
clc

%-------------------------------------------------------------------------
% CUSTOM PARAMETERS (user can change these)

%--------------------
%-- general parameters
fmax = 2000 ; % maximum frequency to be computed
ppw  = 5 ; % points per wavelength at maximum frequency. Choose 3 <= ppw 
Nmodes = [] ; % select total number of modes to be computed. If Nmodes = [], all possible modes are computed
%--------------------

%--------------------
%-- plate parameters
rho = 390 ; % density [kg / m^3]
nux = 0.39 ; % poisson ratio
Ex = 10.4e9 ; % young's mod along x [Pa]
Ey = 0.994e9 ; % young's mod along y [Pa]
Gxy = 0.526e9 ; % shear mod [Pa]
Lx = 0.4 ; % edge lenght x [m]
Ly = 0.6 ; % edge length y [m]
Lz = 3e-3 ; % thickness [m]

%-- elastic constants around the edges
KRmat = [1e13,1e13; %Kx0 Rx0 
    1e13, 1e13; % K0y R0y
    1e13, 1e13; % KxL RxL
    1e13, 1e13] ; % KLy RLy
%--------------------

%--------------------
%-- braces parameters
Nribs = 8 ; % number of braces 

Eb = [10e9,10e9,10e9,10e9,10e9,10e9,10e9,10e9].' ; % youngs moduli
Lzb = [3e-3,3e-3,3e-3,3e-3,3e-3,3e-3,3e-3,3e-3].' ; % thicknesses 
bb = [3e-2,3e-2,3e-2,3e-2,3e-2,3e-2,3e-2,3e-2].' ; % width cross section
rhob = [400,400,400,400,400,400,400,400].' ; % densities

% rib coordinates along x (start and end) AS A FRACTION OF Lx
x_beam_coord = ...
    [0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8] ;

% rib coordinates along y (start and end) AS A FRACTION OF Ly
y_beam_coord = ...
    [0.1,0.2;
    0.2,0.3;
    0.3,0.4;
    0.4,0.5;
    0.5,0.6;
    0.6,0.7;
    0.7,0.8;
    0.8,0.9] ;
%--------------------

%--------------------
%-- static loads and stiffeners parameters

Nlump = 3 ; % number of lumped elements
x_lump_coord = [0.12,0.47,0.91].' ; % x coordinates of lumped elements 
y_lump_coord = [0.25,0.4,0.38].' ; % y coordinates of lumped elements 
KLump  = [1e3,3e5,1.2e2].' ;
MLump  = [0.01,0.01,0.01].' ;
%--------------------

%--------------------
%-- plot parameters parameters
addpath('/Users/micheleducceschi/Documents/MATLAB/DrosteEffect-BrewerMap-3.2.5.0') ; % some cool colormaps 
cmap = brewermap(512, 'PRGn'); % colormap
NN = 1 ; % first mode number to be plotted 
Nplots = 6 ; % select 3,6 or 9. If another number is selected, it is defaulted to 3. Displayed plots are NN + (0:Nplots - 1)
%--------------------

% END CUSTOM PARAMETERS 
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% DERIVED PARAMETERS (user cannot change these)

%-- plate grid spacings and grid points
Evec = [Ex,Ey,Gxy] ;
Lvec = [Lx,Ly,Lz] ;
facx = sqrt(2*pi)*sqrt(sqrt(Evec(1)*Lvec(3)^2/12/rho)) ;
facy = sqrt(2*pi)*sqrt(sqrt(Evec(2)*Lvec(3)^2/12/rho)) ;
hx = facx/ppw/sqrt(fmax) ; % this applies the formula hx*ppw = cx/fmax (\lambda f = c)
hy = facy/ppw/sqrt(fmax) ; % this applies the formula hy*ppw = cy/fmax (\lambda f = c)
Nx   = round(Lvec(1)/hx) ; Ny = round(Lvec(2)/hy) ;
Nvec = [Nx,Ny] ; % grid points
hvec = [Lvec(1)/Nvec(1), Lvec(2)/Nvec(2)] ; % grid spacings 


%-- braces parameters
Ib = bb.*Lzb.^3/12 ; % moments of inertia
Ab = bb.*Lzb ; %  areas
hb = sqrt(sqrt(Eb.*Ib./rhob./Ab))*sqrt(2*pi/fmax)/ppw ; %  grid spacings in terms of ppw

Lb = sqrt((x_beam_coord(:,1)-x_beam_coord(:,2)).^2 + (y_beam_coord(:,1)-y_beam_coord(:,2)).^2) ;
Nb = round(Lb./hb) ; % number of points
hb = Lb./Nb ; % actual grid spacings 

beamParams = [rhob,Eb,Ib,Ab,Nb,hb] ;
beamCoord = [x_beam_coord; y_beam_coord] ;


%-- lumped elements parameters
lumpParams = [KLump,MLump] ;
lumpCoord = [x_lump_coord; y_lump_coord] ;


%-- plot parameters
close all
xax = linspace(0,Lx,Nx+1) ;
yax = linspace(0,Ly,Ny+1) ;
[X,Y] = meshgrid(xax,yax) ;
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% EIGENVALUE PROBLEM
[Om,Q] = freq_domain_sim(rho,Evec,nux,Lvec,hvec,Nvec,KRmat,fmax,Nmodes,Nribs,beamParams,beamCoord,Nlump,lumpParams,lumpCoord) ;
Qplate = Q(1:end-Nlump,:) ;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% PLOT RESULTS
modal_plotter(cmap,Nplots,NN,X,Y,Qplate,Lvec,Nvec,Nribs,Nlump,beamParams,beamCoord,lumpCoord)
%-------------------------------------------------------------------------
